import setuptools
import fompy

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="fti-fompy",
    version=fompy.VERSION,
    author="Arseniy Kononov",
    author_email="a.kononov1@g.nsu.ru",
    description="Routines for the course \"Physical Fundamentals of Microelectronics\"",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kononovarseniy/fompy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
        'scipy'
    ],
    python_requires='>=3.6',
)
