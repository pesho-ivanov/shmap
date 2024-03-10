from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        "sketch",
        ["sketch.pyx"],
        extra_compile_args=[
            "-std=c++20", 
            "-Wno-unreachable-code",
            #"-stdlib=libc++",
            #"-O3", "-ffast-math", "-march=native", "-Xpreprocessor", "-fopenmp"
        ],
        language='c++',
        extra_link_args=[
            #"-stdlib=libc++",
            #"-Xpreprocessor", "-fopenmp"
        ],
        #include_dirs=[numpy.get_include()]
    )
]

setup(
    ext_modules = cythonize(extensions, annotate=True, language_level=3)
)
