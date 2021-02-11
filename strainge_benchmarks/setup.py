from setuptools import setup


setup(
    name='strainge_benchmarks',
    version='1.0',
    py_modules=['strainge_benchmarks', 'longitudinal_truth'],
    scripts=[
        'bin/prepare_strainge_db.py',
        'bin/sample_db.py',
        'bin/benchmark_data_cfg.py',
        'bin/mutate_ref.py',
        'bin/midas_liftover.py',
        'bin/midas_mask_high_density.py',
        'bin/create_contig_map.py',
        'bin/gen_titration_sample.py'
    ],
    url='https://github.com/broadinstitute/strainge_benchmarks',
    author='Bacterial Genomics Lab, Broad Institute',
    license='MIT'
)
