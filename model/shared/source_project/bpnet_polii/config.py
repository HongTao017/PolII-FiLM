from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CONFIG = {
    "run_name": "default",
    "region_lengths": {
        "tss_sequence": 600,
        "first_donor_sequence": 2000,
        "first_acceptor_sequence": 2000,
        "pre_tes_sequence": 1000,
    },
    "filters": 64,
    "conv1_kernel_size": 25,
    "region_dilated_layers": {
        "tss_sequence": 7,
        "first_donor_sequence": 8,
        "first_acceptor_sequence": 8,
        "pre_tes_sequence": 8,
    },
    "regulatory_hidden": 256,
    "regression_hidden": [256, 64],
    "dropout": 0.2,
    "normalization": "groupnorm",
    "groupnorm_groups": 8,
    "use_film": True,
    "batch_size": 64,
    "gradient_accumulation": 1,
    "epochs": 100,
    "patience": 25,
    "lr": 0.0001,
    "wd": 0.00001,
    "film_wd": 0.00001,
    "gradient_clipping": 1.0,
    "num_workers": 4,
    "precision": "bfloat16",
    "checkpoint_metric": "rmse",
    "evaluate_test": True,
    "allow_cross_fold_overlap": False,
    "train_eval_sample_size": 2048,
    "replicate_seeds": [20260915, 20260916, 20260917, 20260918, 20260919],
}
