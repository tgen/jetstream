import os
import io
import yaml


def default():
    return {
        'profile': 'default',
        'backend': 'local',
        'autosave': 0,
        'project_index_dir': 'jetstream',
        'project_temp_dir': 'temp',
        'project_logs_dir': 'logs',
        'project_results_dir': 'results',
        'project_extra_dirs': [],
        'log_level': None,
        'log_format': None,
        'openblas_num_threads': '2',
        'mkl_num_threads': '2',
    }


def profile_path():
    home = os.path.expanduser('~')
    path = os.path.join(home, '.jetstream_profile')
    return os.getenv('JETSTREAM_PROFILE', path)


def load_profile(path=None):
    """Load settings from a yaml file."""
    if path is None:
        path = profile_path()

    settings = default()

    if os.path.isfile(path):
        with open(path, 'r') as fp:
            settings.update(yaml.safe_load(fp))
        settings.update(profile=path)

    return settings


def save_profile(settings, path):
    with open(path, 'w') as fp:
        yaml.dump(settings, indent=4, fp=fp)


def dumps_profile(settings):
    stream = io.StringIO()
    yaml.dump(settings, stream=stream, default_flow_style=False)
    return stream.getvalue()