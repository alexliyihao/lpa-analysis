import os
import pandas as pd
from . import path_settings, utils, settings

class DataManager():
    
    def __init__(self):
        """
        init the DataManagement instance, with a specific version to be used
        """
        # get the root
        self.root = path_settings.ROOT
        
        # load the paths to be defined
        self.paths = utils.DotDict({
            datatype_name: datatype_meta["path"] 
            for datatype_name, datatype_meta 
            in path_settings.DATAMETA.items()
        })
        
        self.filename_formats = utils.DotDict({
            datatype_name: datatype_meta["name_format"] 
            for datatype_name, datatype_meta 
            in path_settings.DATAMETA.items()
        })
        
        # for each datatype, load all the folders under this datatype
        self.versions = utils.DotDict({
            data_type: self._get_subdirectory_list(sub_path = path) 
            for data_type, path 
            in self.paths.items()
        })
        
        # additional metadata
        self._generate_meta()
        
    def _get_subdirectory_list(self, sub_path):
        """
        helper functions for self._generate_meta, get all the folders under a sub path
        args:
            sub_path string, a sub path under self.root
            
        return:
            list[str], a list with all 1 level children under sub_path
        """
        return next(os.walk(os.path.join(self.root, sub_path)))[1]
    
    def _generate_meta(self):
        """
        get the metadata of DataManager, including the versions of all data
        """
        pass
    
    def deploy(self):
        """
        create all the folders for registered data type
        """
        for data_type, path in self.paths.items():
            try:
                os.makedirs(os.path.join(self.root, path), exist_ok = False)
            except FileExistsError:
                print(f"folder for data type {data_type} exists")
                
    def _create_version(self, data_type, version_name):
        """
        base method for all folder building
        """
        path = os.path.join(self.root, self.paths[data_type], version_name)
        try:
            os.makedirs(path, exist_ok = False)
            self.versions[data_type].append(version_name)
            print(f"{data_type[:-1].replace('_', ' ')} version {version_name} created")
        except FileExistsError:
            print(f"folder for data type {data_type} exists")


class DataManager_LPA(DataManager):
    """
    the data manager for LPA pipeline
    """
    def _generate_meta(self):
        """
        a overridden DataManager method generating the data
        """
        # load the metadata folders
        self.coassin_output = path_settings.COASSIN_OUTPUT
        self.paths["analysis_results"] = path_settings.ANALYSIS_RESULT_PATH
        self.versions["analysis_results"] = self._get_subdirectory_list(
            sub_path = path_settings.ANALYSIS_RESULT_PATH)
        
    def create_vcf(self, version_name):
        """
        create folders for new vcf folder with <version_name>, add them to self.versions
        args:
            version_name, string, the name of the version
        return:
            snps_path, genotypes_path, the path of folders saving snps and genotypes
        """
        self._create_version(data_type = "vcfs", version_name = version_name)
    
    def create_dataset(self, version_name):
        """
        create folders for new datasets with <version_name>, add them to self.versions
        args:
            version_name, string, the name of the version
        return:
            snps_path, genotypes_path, the path of folders saving snps and genotypes
        """
        self._create_version(data_type = "snps", version_name = version_name)
        self._create_version(data_type = "genotypes", version_name = version_name)
        
    def create_analysis_result(self, version_name):
        """
        create folders for new analysis_result folder with <version_name>, add them to self.versions
        args:
            version_name, string, the name of the version
        return:
            snps_path, genotypes_path, the path of folders saving snps and genotypes
        """
        self._create_version(data_type = "analysis_results", version_name = version_name)
        
    def _add_prefix(self, parent_list, hispanic_list, i):
        """
        helper function for generate_bam_list, add correct "hispanic" path if necessary
        """
        if i in parent_list:
            return i
        elif i in hispanic_list:
            return os.path.join("hispanic",i)

    def _generate_bam_list(self, parent_path, snp_path):
        """
        the function generate all the relative bam files from a pandas.DataFrame
        Args:
            parent_path: str, the path taking all the coassin output
            snp_path: str, the dataframe to be analyzed
        Return:
            list[str]: the list of all relative bam files to parent path
        """
        df = pd.read_csv(snp_path, index_col = 0)
        parent_list = next(os.walk(parent_path))[1]
        hispanic_list = next(os.walk(os.path.join(parent_path, "hispanic")))[1]
        if "WES_ID" in df.columns:
            return [self._add_prefix(parent_list = parent_list, 
                                    hispanic_list = hispanic_list, 
                                    i = WES_ID+".BQSR.recaled.bam") 
                    for WES_ID in df.WES_ID]
        else:
            return [self._add_prefix(parent_list = parent_list, 
                        hispanic_list = hispanic_list, 
                        i = WES_ID+".BQSR.recaled.bam") 
                    for WES_ID in df.index]

    def path_tool_vcf(self, ethnicity, version_name):
        """
        utility function generate all the related paths about:
        1. coassin output original path
        2. all the coassin output's relative path to coassin output
        3. the name of vcf files
        Args:
            ethnicity: str or int, the ethnicity node
            missing: mode, if missing = True, 20 missing value will be removed from the list
        """
        assert ethnicity in ["eu","af",'hisp',"others", 'complete']
        input_path = self.coassin_output
        bam_list = self._generate_bam_list(
            parent_path = self.coassin_output, 
            snp_path = os.path.join(
                self.root, 
                self.paths["snps"], 
                version_name,
                self.filename_formats["snps"].replace(
                    "*", str(path_settings.CLASS_ALIAS[ethnicity])
                )
            ))
        
        vcf_name = self.filename_formats["vcfs"].replace(
                       "*", str(path_settings.CLASS_ALIAS[ethnicity])
                   )
        
        genotypes_name = self.filename_formats["genotypes"].replace(
                            "*", str(path_settings.CLASS_ALIAS[ethnicity])
                         )
        
        vcf_path = os.path.join(self.root, self.paths["vcfs"],version_name, vcf_name)
        
        genotypes_path = os.path.join(self.root, self.paths["genotypes"], version_name, genotypes_name)
        return utils.DotDict({"input_path":input_path, 
                              "bam_list":bam_list, 
                              "vcf_name":vcf_name, 
                              "vcf_path": vcf_path, 
                              "genotypes_name": genotypes_name, 
                              "genotypes_path": genotypes_path})
    
    def path_tool_csv(self, version_name):
        """
        a wrapper for _generate_paths functions, once called, 
        it will return a dot access dictionary of paths
        """
        snps_path = os.path.join(self.root, self.paths["snps"], version_name)
        genotypes_path = os.path.join(self.root, self.paths["genotypes"], version_name)
        snps_dict = {ethnicity: os.path.join(snps_path, 
                                             self.filename_formats["snps"].replace("*", str(alias))) 
                     for ethnicity, alias 
                     in path_settings.CLASS_ALIAS.items()
                    }
        genotypes_dicts = {ethnicity: os.path.join(snps_path, 
                                                   self.filename_formats["genotypes"].replace("*", str(alias))) 
                           for ethnicity, alias 
                           in path_settings.CLASS_ALIAS.items()
                          }
        return utils.DotDict({"pheno":utils.DotDict(snps_dict),
                              "geno": utils.DotDict(genotypes_dicts)
                             })
        

class DataLoader():
    """
    a data loader for datasets pipeline and vcf pipeline, which can load 1 dataset and one vcf pipeline at maximum
    """
    def __init__(self):
        """
        init function
        """
        self.root = path_settings.ROOT
        self.snps_path = path_settings.SNPS_PATH
        self.genotypes_path = path_settings.GENOTYPES_PATH
        self.vcf_path = path_settings.VCF_PATH
        
    def load_dataset(self, version_name):
        self.snps_dict = {"eu": os.path.join(self.root, self.snps_path, version_name, path_settings.SNPS_NAME_FORMAT.replace("*", "1")),
                          "af": os.path.join(self.root, self.snps_path, version_name, path_settings.SNPS_NAME_FORMAT.replace("*", "2")),
                          "hisp": os.path.join(self.root, self.snps_path, version_name, path_settings.SNPS_NAME_FORMAT.replace("*", "3")), 
                          "complete": os.path.join(self.root, self.snps_path, version_name, path_settings.SNPS_NAME_FORMAT.replace("*", "complete"))
                         }
        self.genotypes_dicts = {"eu": os.path.join(self.root, self.genotypes_path, version_name, path_settings.GENOTYPES_NAME_FORMAT.replace("*", "1")),
                                "af": os.path.join(self.root, self.genotypes_path, version_name, path_settings.GENOTYPES_NAME_FORMAT.replace("*", "2")),
                                "hisp": os.path.join(self.root, self.genotypes_path, version_name, path_settings.GENOTYPES_NAME_FORMAT.replace("*", "3")), 
                                "complete": os.path.join(self.root, self.genotypes_path, version_name, path_settings.GENOTYPES_NAME_FORMAT.replace("*", "complete"))
                               }
    
    def get_dataset(self, ethnicity, n_table):
        assert n_table in [2,3]
        assert ethnicity in ['eu','af','hisp', 'complete']       
        genotype = pd.read_csv(self.genotypes_dicts[ethnicity], index_col = 0)
        snp = pd.read_csv(self.snps_dict[ethnicity], index_col = 0)
        genotype,snp = genotype.align(snp, axis = 0)
        if "AF" in snp.columns:
            snp_x_list = ["GENDER","AGE","AF","HISP"]
        else:
            snp_x_list = ["GENDER","AGE"]
        snp_x = snp[snp_x_list]
        snp_y = snp[settings.Y_LIST]
        if n_table == 2:
            df = pd.concat([genotype, snp_x], axis = 1)
            return df, snp_y
        if n_table == 3:
            return genotype, snp_x, snp_y
        
    def get_snp_ys(self):
        return [pd.read_csv(self.snps_dict[ethnicity], index_col = 0) for ethnicity in ["eu", "af", "hisp", "complete"]]
            
    def _get_vcf_names(self, vcf_path):
        """
        the function get column names from a VCF file, helper function for read_vcf
        Args:
            vcf_path: the path of the files
        Return:
            vcf_names: list[str], the column name of the vcf file
            i: int, the number of lines where the headers are at.
        """
        i = 0
        with open(vcf_path, "rt") as file:
            # keep reading the lines until finding the proper line
            for line in file:
                # the line start with "#CHROM"
                if line.startswith("#CHROM"):
                    vcf_names = [x for x in line.split('\t')]
                    break
                i += 1
        return vcf_names, i

    def _read_vcf(self, path):
        """
        the function reading a simpfies form VCF files into useable csv file
        Args:
            path: str, the path to the file
        Return:
            pandas.DataFrame, the tranposed pandas dataframe can be used in the next step
        """
        # get the column name of the vcf file
        names, i = self._get_vcf_names(path)
        # read the csv file
        vcf = pd.read_csv(path, sep = "\t", header = 0, skiprows = i, low_memory = False)
        # position-ref/alt is the only keys we will use
        vcf["pos-ref/alt"] = vcf.apply(axis = 1, func = lambda row: f"{row.POS}-{row.REF}/{row.ALT}")
        # Drop not-related columns
        vcf.drop(["#CHROM", "ID", "QUAL", "FILTER", "FORMAT", "REF", 'ALT', 'POS'], axis = 1, inplace = True)
        # use pos-ref/alt as the index
        vcf.set_index("pos-ref/alt", inplace = True)
        return vcf.transpose()
    
    def load_vcf(self, version_name):
        assert version_name in settings.VCF_VERSIONS
        path = os.path.join(self.root, self.vcf_path, self.version_name)
        self.vcf_dict = {"eu": os.path.join(self.root, self.snps_path, path_settings.VCFS_NAME_FORMAT.replace("*", "1")),
                         "af": os.path.join(self.root, self.snps_path, path_settings.VCFS_NAME_FORMAT.replace("*", "2")),
                         "hisp": os.path.join(self.root, self.snps_path, path_settings.VCFS_NAME_FORMAT.replace("*", "3")), 
                         "complete": os.path.join(self.root, self.snps_path, path_settings.VCFS_NAME_FORMAT.replace("*", "complete"))
                         }
        
    def get_vcf(self, ethnicity):
        return self._read_vcf(self.vcf_dict[ethnicity])
        


