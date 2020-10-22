import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py_exposure_calculator", # Replace with your own username
    version="0.1",
    author="Mohamed Shaaban",
    author_email="m.shaaban@mail.utoronto.ca",
    description="A simple program for sensitivity and exposure calculations for astronmical missions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/superbit-collaboration/py_exposure_calculator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GPL",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)