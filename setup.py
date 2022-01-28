import setuptools

setuptools.setup(
    packages=['mdvwhole'],
    entry_points = {
        'console_scripts': ['mdvwhole=mdvwhole.whole:main'],
        },
)
