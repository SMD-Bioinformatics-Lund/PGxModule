import os
import pandas as pd
import numpy as np
from datetime import datetime
import jinja2
from weasyprint import HTML


class UtilityFunctions:
    def readfile(self, file_path, sep="\t"):
        try:
            # Attempt to read the file
            df = pd.read_csv(file_path, sep=sep, dtype=str)
            return df
        except FileNotFoundError:
            print(f"Error: File not found at path: {file_path}")
            return None
        except pd.errors.EmptyDataError:
            print(f"Error: File at path {file_path} is empty.")
            return None
        except pd.errors.ParserError:
            print(
                f"Error: Unable to parse the file at path {file_path}. Please check the file format."
            )
            return None

    def risk_duplication(self, x):
        min_dup_var = min(abs(np.array([1 / 4, 1 / 3, 2 / 3, 3 / 4]) - np.array(x)))
        min_diploid = min(abs(np.array([0, 1 / 2, 1]) - np.array(x)))

        if min_dup_var > min_diploid:
            return False
        return True

    def round_df(self, df, digits):
        numeric_columns = df.select_dtypes(include="number").columns
        df[numeric_columns] = df[numeric_columns].round(decimals=digits)
        return df

    def get_haplotypes(self, ID, haplotype_definitions):
        haplo = haplotype_definitions.loc[
            haplotype_definitions["ID"] == ID, "HAPLOTYPE"
        ]
        return "/".join(haplo)


class Report(UtilityFunctions):
    def __init__(self, **kwargs) -> None:
        super().__init__()
        self._check_args(**kwargs)
        self._set_attributes(**kwargs)
        self.targets_depth_df = self.readfile(self.missing_annotated_depth)
        self.haplotype_definitions_df = self.readfile(self.haplotype_definitions)
        self.possible_diplotypes_df = self.readfile(self.possible_diplotypes)
        self.possible_interactions_df = self.readfile(self.possible_interactions)
        self.target_bed_df = self.readfile(self.target_bed)

    def _check_args(self, **kwargs):
        # Check if all required arguments are present
        required_args = [
            "group",
            "read_depth",
            "detected_variants",
            "missing_annotated_depth",
            "haplotype_definitions",
            "possible_diplotypes",
            "possible_interactions",
            "target_bed",
            "padded_baits_depth",
            "target_rsids",
            "annotated_vcf",
            "dbSNP_version",
            "output",
            "report_template",
            "genome_version",
        ]

        missing_args = [arg for arg in required_args if arg not in kwargs]
        if missing_args:
            raise ValueError(f"Missing required arguments: {', '.join(missing_args)}")

    def _set_attributes(self, **kwargs):
        # Set attributes based on kwargs
        for arg_name, arg_value in kwargs.items():
            setattr(self, arg_name, arg_value)

    def get_clinically_relevant_variants(self) -> pd.DataFrame:
        """
        Retrieves clinically relevant variants from the given DataFrame.

        Args:
            variants_df (pd.DataFrame): The DataFrame containing variant information.

        Returns:
            pd.DataFrame: The DataFrame with clinically relevant variants.
        """
        var_df = self.readfile(self.detected_variants)
        if var_df is None:
            return []

        if var_df.empty:
            return []

        # Assuming detected_variants is a DataFrame with a 'GT' column
        var_df["Zygosity"] = var_df["GT"].apply(
            lambda x: "Homo"
            if x in ["0/0", "1/1"]
            else "Hetero"
            if x in ["1/0", "0/1"]
            else None
        )

        var_df["Position"] = var_df.apply(
            lambda row: f"{row['#CHROM']}:{row['POS']}", axis=1
        )
        var_df[["Ref.reads", "Alt.reads"]] = var_df["AD"].str.split(",", expand=True)
        var_df[["Ref.reads", "Alt.reads"]] = var_df[["Ref.reads", "Alt.reads"]].apply(
            pd.to_numeric
        )

        var_df["Haplotype"] = var_df["ID"].apply(
            lambda x: self.get_haplotypes(x, self.haplotype_definitions_df)
        )

        var_df["Variant Frequency"] = var_df["Alt.reads"] / (
            var_df["Alt.reads"] + var_df["Ref.reads"]
        )

        var_df = var_df[var_df["Variant Frequency"] > 0]

        var_df["Possible Duplication"] = var_df["Variant Frequency"].apply(
            lambda x: self.risk_duplication(x)
        )

        var_df["Variant Frequency"] = var_df["Variant Frequency"] * 100
        var_df["Variant Frequency"] = var_df["Variant Frequency"].round(2)

        columns = [
            "GENE",
            "ID",
            "Haplotype",
            "Position",
            "Zygosity",
            "Variant Frequency",
            "Possible Duplication",
        ]

        new_columns = [
            "Gene",
            "rsID",
            "Possible Haplotypes",
            "Position",
            "Zygosity",
            "Variant Frequency",
            "Possible Duplication",
        ]

        variants_present = var_df[columns].copy()
        variants_present = variants_present.rename(
            columns=dict(zip(columns, new_columns))
        )

        return variants_present

    def get_faulty_haplotypes(self, variants_df) -> list:
        """
        Get the faulty haplotypes from the given variants DataFrame.

        Args:
            variants_df (pandas.DataFrame): DataFrame containing the variants.

        Returns:
            list: List of faulty haplotypes.
        """
        if variants_df.empty or variants_df is None:
            return []

        faulty_haplotypes = (
            variants_df.loc[variants_df["Possible Duplication"], "Possible Haplotypes"]
            .str.split("/")
            .explode()
            .unique()
        )

        return faulty_haplotypes

    def get_clinical_recommendations(self, faulty_haplotypes_list) -> pd.DataFrame:
        clin_columns = ["gene", "Haplotype1", "Haplotype2", "Guideline"]
        verbose_columns = [
            "Gene",
            "Haplotype 1",
            "Haplotype 2",
            "Clinical Recommendation",
        ]

        clinical_guidelines_present = self.possible_diplotypes_df[clin_columns].rename(
            columns=dict(zip(clin_columns, verbose_columns))
        )

        faulty_haplotypes_set = set(faulty_haplotypes_list)
        warning_idx = clinical_guidelines_present[
            (clinical_guidelines_present["Haplotype 1"].isin(faulty_haplotypes_set))
            | (clinical_guidelines_present["Haplotype 2"].isin(faulty_haplotypes_set))
        ].index.tolist()

        return clinical_guidelines_present, warning_idx

    def get_interactions_warnings(self, faulty_haplotypes_list) -> list:
        """
        Get the interactions for the given list of faulty haplotypes.

        Args:
            faulty_haplotypes_list (list): A list of faulty haplotypes.

        Returns:
            list: A list of indices indicating the interactions with faulty haplotypes.
        """
        if len(self.possible_interactions_df) != 0:
            interaction_warning_idx = self.possible_interactions_df["haplotypes"].apply(
                lambda x: any(
                    haplotype in faulty_haplotypes_list for haplotype in x.split(",")
                )
            )

            interaction_warning_idx = interaction_warning_idx[
                interaction_warning_idx
            ].index.tolist()
        else:
            interaction_warning_idx = []

        return interaction_warning_idx

    def get_low_depth_targets(self, read_depth_threshold=100) -> pd.DataFrame | None:
        """
        Retrieves the targets with read depth below a specified threshold.

        Args:
            read_depth_threshold (int): The minimum read depth threshold. Targets with read depth below this value will be returned.

        Returns:
            pandas.DataFrame or None: A DataFrame containing the targets with low read depth, or None if no targets meet the criteria.
        """
        columns = ["ID", "Haplotype", "Locus", "Total_Depth"]
        self.targets_depth_df["Haplotype"] = self.targets_depth_df["ID"].apply(
            lambda x: self.get_haplotypes(x, self.haplotype_definitions_df)
        )
        processed_targets_depth_df = self.targets_depth_df.copy()
        processed_targets_depth_df = processed_targets_depth_df[columns]

        colnames = ["rsID", "Haplotype", "Position", "Read Depth"]
        processed_targets_depth_df.columns = colnames
        processed_targets_depth_df["Read Depth"] = pd.to_numeric(
            processed_targets_depth_df["Read Depth"], errors="coerce"
        )

        # Return only the targets with below the depth
        low_depth_targets = processed_targets_depth_df[
            processed_targets_depth_df["Read Depth"] < read_depth_threshold
        ]

        # Display the subset of targets with low depth
        if low_depth_targets.empty:
            return None

        return low_depth_targets

    def df_to_dict(self, df):
        if df is None or df.empty:
            return None
        else:
            return df.to_dict(orient="records")

    def get_targets(self) -> dict | None:
        """
        ID    HAPLOTYPE
        rs1142345    TPMT-3A
        rs1142345    TPMT-3C
        rs1800460    TPMT-3A
        """

        if not self.haplotype_definitions_df.empty:
            target_rsids = (
                self.haplotype_definitions_df[["ID", "HAPLOTYPE"]]
                .drop_duplicates()
                .sort_values("ID")
                .values.tolist()
            )

            target_rsids_dict = {}

            for i in target_rsids:
                if i[0] in target_rsids_dict:
                    target_rsids_dict[i[0]] = f"{target_rsids_dict[i[0]]}/{i[1]}"
                else:
                    target_rsids_dict[i[0]] = i[1]

            return target_rsids_dict
        else:
            return None

    def create_report(self):
        targets = self.get_targets()
        clincial_variants = self.get_clinically_relevant_variants()

        faulty_haplotypes = self.get_faulty_haplotypes(clincial_variants)

        guidelines_present, guidelines_warnings = self.get_clinical_recommendations(
            faulty_haplotypes
        )

        interactions_warnings = self.get_interactions_warnings(faulty_haplotypes)
        low_depth_targets = self.get_low_depth_targets(
            read_depth_threshold=self.read_depth
        )

        # Render the Jinja2 template
        template_loader = jinja2.FileSystemLoader(
            searchpath=os.path.dirname(self.report_template)
        )
        template_env = jinja2.Environment(loader=template_loader)
        template_env.filters["set"] = set
        template_env.filters["list"] = list
        template = template_env.get_template(os.path.basename(self.report_template))

        # Render the template with the data
        rendered_content = template.render(
            group=self.group,
            clincial_variants=self.df_to_dict(clincial_variants),
            faulty_haplotypes=faulty_haplotypes,
            guidelines_present=self.df_to_dict(guidelines_present),
            guidelines_warnings=guidelines_warnings,
            possible_interactions=self.df_to_dict(self.possible_interactions_df),
            interactions_warnings=interactions_warnings,
            low_depth_targets=self.df_to_dict(low_depth_targets),
            targets=targets,
            read_depth=self.read_depth,
            date=datetime.now().strftime("%Y-%m-%d"),
            dbSNP_version=self.dbSNP_version,
            genome_version=self.genome_version,
        )

        # Save the rendered content to the output file
        with open(self.output, "w") as output_file:
            output_file.write(rendered_content)

        # pdf_file_path = "output.pdf"
        # HTML(string=rendered_content).write_pdf(pdf_file_path)

        print(f"Report saved to: {self.output}")
