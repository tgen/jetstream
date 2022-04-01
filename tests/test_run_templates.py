import os
import pytest
import tempfile
import jetstream
import subprocess
from jetstream.cli import main as cli_main

jetstream.settings.clear()
jetstream.settings.read(user=False)
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_TEMPLATES = os.path.join(TESTS_DIR, "templates")

try:
    subprocess.run(
        ["sbatch", "--version"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    is_sbatch_available = True
except FileNotFoundError:
    is_sbatch_available = False

try:
    subprocess.run(
        ["singularity", "--version"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    is_singularity_available = True
except FileNotFoundError:
    is_singularity_available = False

local_tests = {
    "helloworld_1.jst": "single helloworld task",
    "dependencies_1.jst": "dependencies can be declared with before/after",
    "dependencies_2.jst": "dependencies can be declared with input/output",
    "dependencies_3.jst": "tasks with dependencies can be added during run with exec",
    "inheritance_1.jst": "templates can include code from other templates",
    "logging_1.jst": "templates can log to stderr with log global fn",
    "mapping_1.jst": "templates can include properties mapping at the top",
    "retry_1.jst": "tasks can include retry directive that allows tasks to fail and then be run again",
    "stress_1.jst": "runs should not crash due to forking limits",
    "stress_2.jst": "runs should be able to process lots of concurrent tasks",
    "stress_3.jst": "tasks with no command should complete very fast",
    "slurm.jst": "task sbatch_args should override slurm args",
}

slurm_tests = {
    "slurm.jst": "task sbatch_args should override slurm args",
    "helloworld_1.jst": "single helloworld task",
    "dependencies_1.jst": "dependencies can be declared with before/after",
    "dependencies_2.jst": "dependencies can be declared with input/output",
    "dependencies_3.jst": "tasks with dependencies can be added during run with exec",
}

singularity_tests = {
    "singularity_helloworld_1.jst": "single helloworld task using singularity",
}

slurm_singularity_tests = {
    "singularity_helloworld_1.jst": "single helloworld task using singularity and submit to slurm",
}

@pytest.fixture()
def cleandir(tmp_path):
    """test fixture that changes to a clean directory before running"""
    old_dir = os.getcwd()
    os.chdir(tmp_path)
    yield tmp_path
    os.chdir(old_dir)


@pytest.mark.parametrize('template', local_tests.keys(), ids=local_tests.values())
def test_run_template_local(template, cleandir):
    path = os.path.join(TEST_TEMPLATES, template)
    args = ("run", path)
    cli_main(args)


@pytest.mark.skipif(not is_singularity_available, reason="singularity not available")
@pytest.mark.parametrize('template', singularity_tests.keys(), ids=singularity_tests.values())
def test_run_template_slurm(template, cleandir):
    jetstream.settings["backend"] = "slurm"
    path = os.path.join(TEST_TEMPLATES, template)
    args = ("run", path)
    cli_main(args)


@pytest.mark.skipif(not is_sbatch_available, reason="sbatch not available")
@pytest.mark.parametrize('template', slurm_tests.keys(), ids=slurm_tests.values())
def test_run_template_slurm(template, cleandir):
    jetstream.settings["backend"] = "slurm"
    path = os.path.join(TEST_TEMPLATES, template)
    args = ("run", path)
    cli_main(args)


@pytest.mark.skipif(not is_sbatch_available and not is_singularity_available, reason="sbatch not available")
@pytest.mark.parametrize('template', slurm_singularity_tests.keys(), ids=slurm_singularity_tests.values())
def test_run_template_slurm(template, cleandir):
    jetstream.settings["backend"] = "slurm"
    path = os.path.join(TEST_TEMPLATES, template)
    args = ("run", path)
    cli_main(args)
