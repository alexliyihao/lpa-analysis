import pandas as pd
import os
from . import settings

class Validators():
    
    def __init__(self, verbose):
        self.verbose = verbose
    
    def _get_folders(self, csv_path):
        """
        get a dictionary of folders need for a specific csv output
        """
        df = pd.read_csv(csv_path, index_col = 0)
        WES_ID = dict(df.WES_ID)
        return {key:f"{value}.BQSR.recaled.bam" for key, value in WES_ID.items()}
    
    def validate(self, csv_file):
        bam_dict = self._get_folders(csv_file)
        bad_dict = self._list_invalid_input(bam_dict)
        if len(bad_dict) == 0:
            if self.verbose:
                print(f"Pass: All the related files to {os.path.basename(csv_file)} is validated")
            return []
        else:
            if self.verbose:
                print(f"Failed: the following {len(bad_dict)} items need a further check")
                for key, values in bad_dict.items():
                    print(f"row index {key}, file {values}")
            return bad_dict
    
    def _validate_input(self, path):
        pass
    
    def _list_invalid_input(self, bam_dict):
        """
        run _validate_input on all the returnings
        """
        return {key: value for key, value in bam_dict.items() if not self._validate_input(value)}
    
    
class SourceValidator(Validators):
    
    def _validate_input(self, path):
        """
        give a path, check if it's a exist path, and if the Coassin pipeline output is there
        """
        return (os.path.isfile(os.path.join(settings.WHICAP_SOURCE, path))
                or 
                os.path.isfile(os.path.join(settings.WHICAP_SOURCE, "hispanic", path)))
    
class OutputValidator(Validators):

    def _validate_input(self, path):
        """
        give a path, check if it's a exist path, and if the Coassin pipeline output is there
        """
        return ((os.path.isdir(os.path.join(settings.COASSIN_OUTPUT, path)) 
                and os.path.isfile(os.path.join(settings.COASSIN_OUTPUT, path, "variantsAnnotate/variantsAnnotate.txt")))
                or 
                (os.path.isdir(os.path.join(settings.COASSIN_OUTPUT, "hispanic", path)) 
                and os.path.isfile(os.path.join(settings.COASSIN_OUTPUT, "hispanic", path, "variantsAnnotate/variantsAnnotate.txt")))
               )
    
class EnsembleValidator():
    def __init__(self):
        self.SV = SourceValidator(verbose = False)
        self.OV = OutputValidator(verbose = False)
    
    def validate(self,path):
        source_bad = self.SV.validate(path)
        output_bad = self.OV.validate(path)
        if (len(source_bad) == 0 and len(output_bad) == 0):
            print(f"Pass: All the related files to {os.path.basename(path)} is validated")
            return (0,0)
        elif (len(source_bad) == 0 and len(output_bad) != 0):
            print(f"Source file Passed")
            print(f"Output File Failed: the following {len(output_bad)} folders need a further check")
            for key, values in output_bad.items():
                print(f"row index {key}, folder {values}")
            return (0,output_bad)
        else:
            output_bad = {key:value for(key, value) in output_bad.items() if (key,value) not in source_bad.items()}
            print(f"Source file Failed: the following {len(source_bad)} files need a further check")
            for key, values in source_bad.items():
                print(f"row index {key}, file {values}")
            if len(output_bad) == 0:
                print("Output files Passed")
                return (source_bad, 0)
            else:
                print(f"Output File Failed: the following {len(output_bad)} folders need a further check")
                for key, values in output_bad.items():
                    print(f"row index {key}, folder {values}")
                return (source_bad, output_bad)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", help="The csv path of folder storing phenotypes")
    
    EV = EnsembleValidator()
    source_bad, output_bad = EV.validate(path = args.path)