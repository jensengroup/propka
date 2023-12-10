import propka.input as m
import zipfile


def test_open_file_for_reading(tmp_path):
    path = tmp_path / "tmp.txt"
    path.write_text("One\nTwo\nThree\n")
    # str
    with m.open_file_for_reading(str(path)) as outer:
        assert outer.read() == "One\nTwo\nThree\n"
    assert outer.closed
    # Path
    with m.open_file_for_reading(path) as outer:
        # TextIO
        with m.open_file_for_reading(outer) as inner:
            assert inner.readline() == "One\n"
        assert not outer.closed
        assert outer.readline() == "Two\n"
    assert outer.closed


def test_open_file_for_reading__zipfile(tmp_path):
    zippath = tmp_path / "tmp.zip"
    arcname = "foo/bar.txt"
    with zipfile.ZipFile(zippath, "w") as ziphandle:
        ziphandle.writestr(arcname, "One\nTwo\nThree\n")
    with m.open_file_for_reading(zippath / arcname) as outer:
        assert outer.readline() == "One\n"
    assert outer.closed
