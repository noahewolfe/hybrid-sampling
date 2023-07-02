from setuptools import setup


def get_requirements():
    with open("requirements.txt", "r") as ff:
        requirements = ff.readlines()
    return requirements


setup(
    name="hybrid_sampling",
    version="0.1.0",
    install_requires=get_requirements(),
    packages=["hybrid_sampling"],
    entry_points={
        "console_scripts": [
            "bilby_pipe_hybrid=hybrid_sampling.hybrid_post:main",
            "bilby_pipe_build_hybrid=hybrid_sampling.build_hybrid:main",
        ]
    },
)
