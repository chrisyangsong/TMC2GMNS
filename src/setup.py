import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tmc2gmns", 
    version="0.1.4",
    author="Dr. Xuesong Zhou, Xianbiao Hu, Jiawei Lu, Chris Yang Song",
    author_email="xzhou74@asu.edu, xbhu@mst.edu, jiaweil9@asu.edu, 845107945@qq.com",
    description="An open-source python package converting TMC data into GMNS format.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[                                 
    # 'os',
    # 'datetime',
    'numpy>=1.19.2',
    'pandas>=1.1.3',
    'MapMatching4GMNS>=0.2.11'
    ],
    url="https://github.com/vegampire/TMC2GMNS",
    packages=[''],
    package_dir={'': 'src'},
    package_data={'': ['src/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)

