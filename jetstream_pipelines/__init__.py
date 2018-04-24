import pkg_resources
__version__ = pkg_resources.get_distribution("jetstream").version

from jetstream_pipelines.jinja import env, package_loader
