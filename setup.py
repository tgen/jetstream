import os
from setuptools import setup, find_packages

__package__ = 'jetstream'
src_dir = os.path.join(os.path.dirname(__file__))
pkg_init = os.path.join(src_dir, __package__, '__init__.py')
readme_path = os.path.join(src_dir, 'README.md')

 
with open(pkg_init) as fp:
    for line in fp:
        if line.startswith('__version__'):
            __version__ = eval(line.split('=', 1)[1])
            break


with open(readme_path) as fp:
    __readme__ = fp.read()



setup(
    name=__package__,
    version=__version__,
    author='Ryan Richholt',
    author_email='ryan@tgen.org',
    url='https://github.com/tgen/jetstream',
    description='NGS analysis pipeline at TGen.',
    long_description=__readme__,
    long_description_content_type='text/markdown',
    keywords='ngs pipeline automation',
    packages=find_packages(exclude=('test',)),
    include_package_data=True,
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
    ],
    install_requires=[
        'networkx',
        'pyyaml',
        'ulid-py',
        'jinja2',
        'filelock',
        'confuse'
    ],
    package_data={
        'jetstream': ['config_default.yaml']
    },
    entry_points={
        'console_scripts': [
            'jetstream=jetstream.cli:main',
        ],
    }
)
