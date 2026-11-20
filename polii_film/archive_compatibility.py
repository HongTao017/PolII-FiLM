"""Restore immutable data releases into the current English source layout.

Transfer checksums describe the original assets. Revision records describe text
metadata after localization; array and checkpoint bytes are never rewritten.
"""
from pathlib import Path
import base64
import hashlib
import json
import re
from functools import lru_cache

ROOT = Path(__file__).resolve().parents[1]
HASH = re.compile(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])')
TEXT_EXTENSIONS = {'.json', '.jsonl', '.md', '.txt', '.tsv', '.py', '.sh', '.r', '.mjs'}
IDENTITY_FIELDS = {'protocol_sha', 'archive_protocol_sha256'}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=4)
def configuration(root=ROOT):
    return json.loads((Path(root) / 'configs/archive_compatibility.json').read_text(encoding='utf-8'))


@lru_cache(maxsize=4)
def replacements(root=ROOT):
    # Legacy UTF-8 tokens are encoded so the public configuration is ASCII.
    return [(base64.b64decode(row['legacy_utf8_base64']).decode('utf-8'), row['replacement'])
            for row in configuration(root)['text_replacements']]


def normalize_text(text, root=ROOT):
    for old, new in replacements(root):
        text = text.replace(old, new)
    return text


def normalize_path(name, root=ROOT):
    # Only path tokens apply here; prose translations cannot rename data files.
    for row in configuration(root)['path_replacements']:
        name = name.replace(base64.b64decode(row['legacy_utf8_base64']).decode('utf-8'), row['replacement'])
    return name


def state_path(root=ROOT):
    return Path(root) / '.cache/archive_compatibility_state.json'


def load_state(root=ROOT):
    path = state_path(root)
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else {'resources': {}, 'protocol_revisions': {}}


def preserve_metadata(relative, content, root=ROOT):
    backup = Path(root) / '.cache/archive_original_metadata' / relative
    backup.parent.mkdir(parents=True, exist_ok=True)
    if not backup.exists():
        backup.write_bytes(content)


def align_metadata(root=ROOT):
    """Align only downloaded text metadata and record each derived revision."""
    root = Path(root)
    state = load_state(root)
    source = json.loads((root / 'configs/source_revision_map.json').read_text(encoding='utf-8'))
    candidates = {}
    for row in source['files']:
        candidates.setdefault(row['original_sha256'], set()).add(row['current_sha256'])
    baseline = {old: next(iter(values)) for old, values in candidates.items() if len(values) == 1}
    paths = [root / rel for rel in state['resources'] if (root / rel).is_file()]
    originals = {}
    for path in paths:
        backup = root / '.cache/archive_original_metadata' / path.relative_to(root)
        originals[path] = (backup if backup.exists() else path).read_bytes()
    # A dependency can refer to a metadata file that itself refers to a source
    # file. Resolve these chains without changing frozen result identities.
    for _ in range(20):
        revisions = dict(baseline)
        for path in paths:
            record = state['resources'][path.relative_to(root).as_posix()]
            revisions[record['original_sha256']] = sha256(path)
        changed = False
        for path in paths:
            content = path.read_bytes()
            try:
                text = content.decode('utf-8-sig')
            except UnicodeError:
                continue
            previous = text
            # Start from the immutable original on every dependency pass. This
            # prevents stale intermediate hashes in multi-level manifests.
            text = normalize_text(originals[path].decode('utf-8-sig'), root)
            def convert(value, key=''):
                if key in IDENTITY_FIELDS:
                    return value
                if isinstance(value, dict):
                    return {k: convert(v, k) for k, v in value.items()}
                if isinstance(value, list):
                    return [convert(v) for v in value]
                if isinstance(value, str):
                    return HASH.sub(lambda m: revisions.get(m.group(), m.group()), value)
                return value
            if path.suffix.lower() == '.json':
                obj = json.loads(text)
                text = json.dumps(convert(obj), ensure_ascii=True, indent=2) + '\n'
            else:
                text = HASH.sub(lambda m: revisions.get(m.group(), m.group()), text)
            if text != previous:
                preserve_metadata(path.relative_to(root), content, root)
                path.write_text(text, encoding='utf-8', newline='\n')
                changed = True
        if not changed:
            break
    else:
        raise ValueError('Metadata revision dependencies did not converge')
    for path in paths:
        state['resources'][path.relative_to(root).as_posix()]['current_sha256'] = sha256(path)
    state_path(root).parent.mkdir(parents=True, exist_ok=True)
    state_path(root).write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')


def record_resource(path, original_bytes, root=ROOT):
    """Localize a small downloaded text member and retain its original bytes."""
    path = Path(path)
    if path.suffix.lower() not in TEXT_EXTENSIONS or len(original_bytes) > 16 * 1024 * 1024:
        return
    try:
        text = original_bytes.decode('utf-8-sig')
    except UnicodeError:
        return
    relative = path.relative_to(root).as_posix()
    text = normalize_text(text, root)
    if path.suffix.lower() == '.json':
        # Escape remaining localized historical labels without changing values.
        obj = json.loads(text)
        text = json.dumps(obj, ensure_ascii=True, indent=2) + '\n'
    if text.encode('utf-8') != original_bytes:
        preserve_metadata(relative, original_bytes, root)
        path.write_text(text, encoding='utf-8', newline='\n')
    state = load_state(root)
    state['resources'][relative] = {'original_sha256': hashlib.sha256(original_bytes).hexdigest(),
                                    'current_sha256': sha256(path)}
    state_path(root).parent.mkdir(parents=True, exist_ok=True)
    state_path(root).write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')


def protocol_fingerprint(path):
    """Return the identity embedded in frozen arrays for a registered revision.

    A newly prepared protocol without an archive identity uses its actual SHA.
    A changed archive protocol must match a registered current content digest.
    """
    path = Path(path)
    actual = sha256(path)
    obj = json.loads(path.read_text(encoding='utf-8'))
    original = obj.get('archive_protocol_sha256')
    if original is None:
        return actual
    known = json.loads((ROOT / 'configs/source_revision_map.json').read_text(encoding='utf-8'))['protocol_revisions']
    known.update(load_state()['protocol_revisions'])
    if known.get(actual) != original:
        raise ValueError('Archive protocol differs from its registered metadata revision')
    return original
