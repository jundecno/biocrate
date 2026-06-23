def test_top_level_import_is_lightweight():
    import biocrate

    assert isinstance(biocrate.__version__, str)


def test_metrics_submodule_imports():
    from biocrate.metrics import accuracy_score, evaluate_regression

    assert accuracy_score([1, 0, 1], [1, 1, 1]) == 2 / 3
    assert set(evaluate_regression([1, 2, 3], [1, 2, 4])) == {
        "mae",
        "mse",
        "rmse",
        "r2",
        "mape",
        "pearsonr",
        "spearmanr",
    }
