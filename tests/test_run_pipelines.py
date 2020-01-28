import os
import tempfile
import pytest
import jetstream
import jetstream.cli

jetstream.settings.clear()
jetstream.settings.read(user=False)
tests_dir = os.path.dirname(os.path.abspath(__file__))
tests_pipelines_dir = os.path.join(tests_dir, 'pipelines')
jetstream.settings['pipelines']['searchpath'] = tests_pipelines_dir


@pytest.fixture()
def cleandir():
    """test fixture that changes to a clean directory before running"""
    old_dir = os.getcwd()
    tempdir = tempfile.TemporaryDirectory()
    os.chdir(tempdir.name)
    yield tempdir
    tempdir.cleanup()
    os.chdir(old_dir)


def main(*args):
    """run the jetstream cli main with given args"""
    jetstream.cli.main(args)


def test_foopipe_1(cleandir):
    """the most basic pipeline requires pipeline.yaml and main template"""
    main('pipelines', 'foopipe_1')


def test_foopipe_2(cleandir):
    """other pipeline features include bin directory"""
    main('pipelines', 'foopipe_2')


def test_foopipe_3(cleandir):
    """pipelines nested in other pipelines can be discovered as well"""
    # barpipe_1 is nested in foopipe_3
    main('pipelines', 'barpipe_1')


def test_foopipe_4(cleandir):
    """pipeline variables are exported and allow files to be accessed"""
    main('pipelines', 'foopipe_4')


def test_resetpipe_1(cleandir):
    """Test resetpipe by running twice without saving progress"""
    with pytest.raises(SystemExit) as e:
        # We expect this to fail the first time
        main('pipelines', 'resetpipe')

    # then pass the second time
    main('pipelines', 'resetpipe')


def test_resetpipe_2(cleandir):
    main('init')

    """Test resetpipe by running twice with saving progress"""
    with pytest.raises(SystemExit) as e:
        # We expect this to fail the first time
        main('pipelines', 'resetpipe')

    # then pass the second time
    main('pipelines', 'resetpipe')
