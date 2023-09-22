from setuptools import setup

setup(
    name="casadi-bspline",
    version=1.0,
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "matplotlib",
        "tk",
        "scipy",
        "wheel",
        "casadi",
    ],
    dependency_links = [],
    packages = [
        'casadi_bspline'
    ],
    package_dir = {
        'casadi_bspline': 'src'
    },
    entry_points = {},
    author="Maksim Surov",
    author_email="surov.m.o@gmail.com",
    description="casadi symbolic splines",
)
