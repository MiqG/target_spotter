import pandas as pd
import defaults
from SplicingDependency import SplicingDependency, PredictFromFile
import utils
import streamlit as st
import tempfile
import os
from datetime import datetime
import shutil
import base64

LOGO_FILE = defaults.LOGO_FILE
EXAMPLE_FILES = defaults.EXAMPLE_FILES
N_JOBS = 5

##### FUNCTIONS #####
def make_download_link(label, source_filename, downloaded_filename):
    with open(source_filename, "rb") as f:
        bytes = f.read()
        b64 = base64.b64encode(bytes).decode()
        download_link = f"<a href='data:file/zip;base64,{b64}' download='{downloaded_filename}'> {label} </a>"

    return download_link


def make_sidebar():

    st.sidebar.title("Target Spotter")
    st.sidebar.image(LOGO_FILE, width=200)
    st.sidebar.markdown(
        """
    Welcome to `target_spotter` App.

    Here, you will be able to transform your splicing and gene expression profiles into splicing dependencies that uncover the weak spots of your cancer samples.
    """
    )

    st.sidebar.markdown(
        """
    ## Contact
This project has been fully developed at the [Centre for Genomic Regulation (CRG)](https://www.crg.eu/) within the group of [Design of Biological Systems](https://www.crg.eu/en/luis_serrano).

Please, report any issues that you experience through this repository's ["Issues"](https://github.com/CRG-CNAG/target_spotter/issues) or email:
- [Miquel Anglada-Girotto](mailto:miquel.anglada@crg.eu)
- [Luis Serrano](mailto:luis.serrano@crg.eu)
    """
    )

    st.sidebar.markdown(
        """
    ## License

`target_spotter` is distributed under a BSD 3-Clause License (see [LICENSE](https://github.com/CRG-CNAG/target_spotter/blob/main/LICENSE)).

    """
    )

    st.sidebar.markdown(
        """
    ## References
- *Himberg, J., & Hyvarinen, A.* "Icasso: software for investigating the reliability of ICA estimates by clustering and visualization". IEEE XIII Workshop on Neural Networks for Signal Processing (2003). DOI: https://doi.org/10.1109/NNSP.2003.1318025
    """
    )

    st.sidebar.markdown("## FAQ")
    with st.sidebar.expander("What are splicing dependencies?", expanded=False):
        st.markdown(
            "Splicing dependencies are calculated from mRNA and splicing profiles of cancer samples"
        )
    with st.sidebar.expander("How should I format my inputs?", expanded=False):
        download_link = "https://drive.google.com/uc?export=download&id=1kHOnt5aDKcxHSuq6bZji9FrH0aQInooi"
        st.markdown(
            """
        To estimate splicing dependencies, `target_spotter` requires two input tables:
        - splicing: PSI profiles of your samples. The first column must contain the splicing events identifiers of [`VastDB`](https://vastdb.crg.eu/wiki/Main_Page) for the Hs2 event annotation (e.g. HsaEX0066564).
        - gene expression: mRNA levels profiles of your samples. The first column must contain the ENSEMBL gene identifiers (e.g. ENSG00000141510)
        
        In case you'd like to see how sample data tables look like, click [here](%s) to download them.
        """
            % download_link
        )


def select_options():
    # inputs
    with st.expander("Upload your own data", expanded=False):
        st.info(
            """
            - Upload your .tsv / .tsv.gz file here. Maximum file size is 200 Mb.
            - Each row corresponds to a molecular feature, each column to a sample.
            - The first column contains the row identifiers.
        """
        )
        uploaded_splicing_file = st.file_uploader(
            "Upload your splicing dataset below", type=["tsv", "tsv.gz"]
        )
        uploaded_genexpr_file = st.file_uploader(
            "Upload your gene expression dataset below", type=["tsv", "tsv.gz"]
        )

        st.markdown(
            """**Note:** 
            By uploading a file, you agree to our
            [Apache License](https://github.com/CRG-CNAG/target_spotter/LICENSE).
            Data that is uploaded via the file uploader will not be saved by us;
            it is only stored temporarily in RAM to perform the calculations.
            """
        )

    # parameters
    with st.expander("Set parameters", expanded=True):
        genexpr_units = st.selectbox(
            "Gene expression units:",
            ["TPM", "Counts"],
        )
        if genexpr_units == "TPM":
            normalize_counts = False
        elif genexpr_units == "Counts":
            normalize_counts = True

    # sample dataset
    with st.expander("Use a sample dataset", expanded=False):
        sample_dataset = st.selectbox("", ["None", "CCLE"])

    if sample_dataset != "None":
        uploaded_splicing_file = EXAMPLE_FILES[sample_dataset]["splicing"]
        uploaded_genexpr_file = EXAMPLE_FILES[sample_dataset]["genexpr"]

    # prepare outputs
    files = {"splicing": uploaded_splicing_file, "genexpr": uploaded_genexpr_file}

    params = {"normalize_counts": normalize_counts}

    return files, params


def upload_data(files, params):
    dirpath = os.path.join(tempfile.mkdtemp(), "splicing_dependency")
    predictor = PredictFromFile(
        splicing_file=files["splicing"],
        genexpr_file=files["genexpr"],
        output_dir=dirpath,
        normalize_counts=params["normalize_counts"],
        n_jobs=N_JOBS,
    )
    predictor.load_data()

    st.info("Data loaded successfully!")

    return predictor


def estimate_dependencies(predictor):
    st.info("Estimating splicing dependencies...")
    estimator = SplicingDependency(
        normalize_counts=predictor.normalize_counts,
        n_iterations=predictor.n_iterations,
        n_jobs=predictor.n_jobs,
    )
    _ = estimator.predict(predictor.splicing_, predictor.genexpr_,)
    st.info("Finished!")
    return estimator


def save_outputs(predictor, estimator):
    # save outputs
    predictor.save(estimator)

    # zip outputs
    output_dir = predictor.output_dir
    shutil.make_archive(output_dir, "zip", output_dir)

    # create link to download outputs
    to_download = output_dir + ".zip"
    now = datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = "splicing_dependency-%s.zip" % now
    download_link = make_download_link(
        "Download Splicing Dependency",
        source_filename=to_download,
        downloaded_filename=filename,
    )
    return download_link


def main():
    st.set_page_config(
        page_title="Target Spotter",
        page_icon=LOGO_FILE,
        layout="centered",
        initial_sidebar_state="auto",
    )

    make_sidebar()

    st.title("Target Spotter")
    st.markdown("Systematic prioritization of splicing targets to treat cancer.")

    left_column, right_column = st.columns(2)
    with left_column:
        st.markdown("## Start your analysis")
        files, params = select_options()
        if st.button("Run Analysis") & (
            (files["splicing"] is not None) | (files["genexpr"] is not None)
        ):
            predictor = upload_data(files, params)
            estimator = estimate_dependencies(predictor)
            download_link = save_outputs(predictor, estimator)
        else:
            download_link = ""
            st.warning("Please, upload your own data or select a sample dataset")

    with right_column:
        st.markdown(download_link, unsafe_allow_html=True)


##### SCRIPT #####
if __name__ == "__main__":
    main()
