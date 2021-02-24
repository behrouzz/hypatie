import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gravitational",
    version="0.0.1",
    author="Behrouz Safari",
    author_email="behrouz.safari@gmail.com",
    description="A python package for querying NASA's Horizons API",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/behrouzz/hypatie",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["hypatie"],
    include_package_data=True,
    install_requires=["numpy"],
    python_requires='>=3.6',
)
