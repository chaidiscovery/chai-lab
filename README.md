# Chai-1

Chai-1 is a multi-modal foundation model for molecular structure prediction that performs at the state-of-the-art across a variety of benchmarks. Chai-1 enables unified prediction of proteins, small molecules, DNA, RNA, glycosylations, and more.

<p align="center">
    <img src='https://github.com/chaidiscovery/chai-lab/blob/main/assets/performance_barplot.png' >
</p>

For more information on the model's performance and capabilities, see our [technical report](https://www.biorxiv.org/content/10.1101/2024.10.10.615955).

## Installation

```shell
# version on pypi:
pip install chai_lab==0.4.1

# newest available version (updates daily to test features that weren't released yet):
pip install git+https://github.com/chaidiscovery/chai-lab.git
```

This Python package requires Linux, and a GPU with CUDA and bfloat16 support. We recommend using an A100 80GB or H100 80GB chip, but A10s and A30s should work for smaller complexes. Users have also reported success with consumer-grade RTX 4090.

## Running the model

The model accepts inputs in the FASTA file format, and allows you to specify the number of trunk recycles and diffusion timesteps via the `chai_lab.chai1.run_inference` function. By default, the model generates five sample predictions, and uses embeddings without MSAs or templates.

The following script demonstrates how to provide inputs to the model, and obtain a list of PDB files for downstream analysis:

```shell
python examples/predict_structure.py
```

To get the best performance, we recommend running the model with MSAs. The following script demonstrates how to provide MSAs to the model by calling out to an MSA server:

```shell
python examples/msas/predict_with_msas.py
```

For further instructions, see `"How can MSAs be provided to Chai-1?"` below.

<details>
<summary>Where are downloaded weights stored?</summary>
<p markdown="1">
By default, weights are automatically downloaded and stored in <package_root>/downloads (usually that's within site-packages).
In cases where you want to control the download location (e.g. on a mounted drive in Docker), you can use the CHAI_DOWNLOADS_DIR envvar to control the download location. For example:

```bash
CHAI_DOWNLOADS_DIR=/tmp/downloads python ./examples/predict_structure.py 
```
</p>
</details>

<details>
<summary>How can MSAs be provided to Chai-1?</summary>
<p markdown="1">

Chai-1 supports MSAs provided as an `aligned.pqt` file. This file format is similar to an `a3m` file, but has additional columns that provide metadata like the source database and sequence pairing keys. We provide code to convert `a3m` files to `aligned.pqt` files. For more information on how to provide MSAs to Chai-1, see [this documentation](examples/msas/README.md).

</p>
</details>

<details>
<summary>How can I customize the inputs to the model further?</summary>
<p markdown="1">

For more advanced use cases, we also expose the `chai_lab.chai1.run_folding_on_context`, which allows users to construct an `AllAtomFeatureContext` manually. This allows users to specify their own templates, MSAs, embeddings, and constraints. We currently provide an example of how to construct an embeddings context as well as an MSA context, and will be releasing helper methods to build template contexts soon.

</p>
</details>

## ⚡ Try it online

We provide a [web server](https://lab.chaidiscovery.com) so you can test the Chai-1 model right from your browser, without any setup.

<p align="center">
    <img src='assets/chailab_online_screenshot.png' height=400 >
</p>

## Using experimental restraints
Chai-1 uniquely offers the ability to fold complexes with user-specified "restraints" as inputs. These restraints specify inter-chain contacts at various resolutions that are used to guide Chai-1 in folding the complex. See [restraints documentation](examples/restraints/README.md) for details.

<p align="center">
    <img src='assets/chailab_restraints_screenshot.png' height=400 >
</p>

## 💬 Feedback

Found a 🐞? Please report it in GitHub [issues](https://github.com/chaidiscovery/chai-lab/issues).

We welcome community testing and feedback. To share observations about the model's performance, please reach via [GitHub discussions](https://github.com/chaidiscovery/chai-lab/discussions), or [via email](mailto:feedback@chaidiscovery.com).

## 🛠️ Development

We use [devcontainers](https://code.visualstudio.com/docs/devcontainers/containers) in development, which helps us ensure we work in identical environments. We recommend working inside a devcontainer if you want to make a contribution to this repository.

Devcontainers work on local Linux setup, and on remote machines over an SSH connection.

## Status

Since this is an initial release, we expect to make some breaking changes to the API and are not guaranteeing backwards compatibility. We recommend pinning the current version in your requirements, i.e.:

```
chai_lab==0.4.1
```

## Citations

If you find Chai-1 useful in your research or use any structures produced by the model, we ask that you cite our technical report:

```
@article{Chai-1-Technical-Report,
	title        = {Chai-1: Decoding the molecular interactions of life},
	author       = {{Chai Discovery}},
	year         = 2024,
	journal      = {bioRxiv},
	publisher    = {Cold Spring Harbor Laboratory},
	doi          = {10.1101/2024.10.10.615955},
	url          = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.10.615955},
	elocation-id = {2024.10.10.615955},
	eprint       = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.10.615955.full.pdf}
}
```

## Licence 

Chai-1 is released under an Apache 2.0 License, which means it can be used for both academic
and commerical purposes, including for drug discovery.

See [LICENSE.md](LICENSE.md).

To discuss partnership and access to new internal capabilities, reach us [via email](mailto:partnerships@chaidiscovery.com).
