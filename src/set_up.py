import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tmc2gmns", 
    version="0.1.0",
    author="Dr. Xuesong Zhou, Xianbiao Hu, Jiawei Lu, Chris Yang Song",
    author_email="xzhou74@asu.edu, xbhu@mst.edu, jiaweil9@asu.edu, songyang0714@gmail.com",
    description="An open-source python package converting TMC data into GMNS format.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vegampire/TMC2GMNS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
