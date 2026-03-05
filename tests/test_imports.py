import lineageresolver as lr


def test_package_import_smoke() -> None:
    assert callable(lr.infer)
    assert callable(lr.estimate_ambient_profile)
