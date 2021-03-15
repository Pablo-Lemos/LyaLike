# example cobaya-compliant SO likelihood package;
# adapted from github.com/cobayasampler/example_external_likelihood

from setuptools import setup

setup(
    name="lyalike",
    version="0.0",
    description="Lyman alpha Likelihoods & Theories",
    zip_safe=False,
    packages=["lyalike"],
    install_requires=[
        "cobaya @ git+https://github.com/cobayasampler/cobaya",
        ],
    #test_suite="soliket.tests",
)
