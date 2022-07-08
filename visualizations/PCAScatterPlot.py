import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


class PCAScatterPlot():
    """A tool class running scatterplot from EIGENSOFT smartPCA result

    For dependency easiness, this class didn't uses ..EIGEN.EIGENToolkit.EIGENToolkit
    instance but copied its read_PCA_result method, feel free to use it individually

    Initialize:
        psp = PCAScatterPlot(output_path = "some_path",
                             display = False,
                             ethnicities = ["eu", "af", "hisp"],
                             **kwargs)

    Plot:
        psp.PCAplot(csv_path, x, y, title,
                    output_path = None, display = False, **kwargs)

        psp.PCAplot_one_vs_others(csv_path, x, y, title, one,
                    output_path = None, display = False, **kwargs)

        output_path: str, the saving folder(not the path of file) of the figure,
                     if provided, it will cover the previously provided path
                     JUST IN THIS RUN
        title: str, the title of the figure, and also the saving file name
               of the figure
        display: Optional[bool], if provided, it will cover the previously
                provided display setting JUST IN THIS RUN
        kwargs: any argument to be passed to sns.scatterplot method,
                for plotting convenience, these setting will PERMANANTLY cover
                the previous setting

    """

    def __init__(self,
                 ethnicities = ["eu", "af", "hisp"],
                 output_path = None,
                 display = False,
                 **kwargs):
        """
        Args:
            output_path: Optional[string], the path saving the image if provided,
                        otherwise the figure will not be saved
            ethnicities: Optional[list[str]], the ethnicities to be take care of,
                        only useful when using self.PCAplot_one_vs_others method
            display: Optional[bool], default False: if True, the image will be
                     displayed once plotted.
            kwargs: any argument to be passed to sns.scatterplot method
        """
        self._ethnicities = ethnicities
        self.output_path = output_path
        self._display = display
        self._sns_setting = kwargs

    def _plot(self, data, x, y, title, output_path = None, display = False):
        """function running the actual visualization from sns.scatterplot

        Args:
            data: pd.DataFrame, the DataFrame to be plotted, which should have
                  a "ethnicity" column for the coloring
        """

        # plot the scatterplot
        sns.scatterplot(data = data,
                        x = x,
                        y = y,
                        hue = "ethnicity",
                        **self._sns_setting).set(title = title)
        # if the output_path is not given, use the original output path
        if output_path is None:
            output_path = self.output_path
        # if output_path is given for any one time, save it
        if not output_path is None:
            saving_path = os.path.join(output_path, f"{title}.png")
            plt.savefig(saving_path)
            print(f"figure saved in {saving_path}")
        # similarly, if the display is not given, use the original display setting
        if display is None:
            display = self._display
        # if the display is True, plot it
        if display == True:
            plt.show()

    def _read(self, path):
        """funtion reading the EIGENSOFT smartPCA results"""
        eigenstrat_result = pd.read_csv(path,
                                 skiprows = 1,
                                 sep = "\s+",
                                 header = None)
        eigenstrat_result = eigenstrat_result.rename(columns = {0: "ID",
                                                                1: "PC1",
                                                                2: "PC2",
                                                                3: "PC3",
                                                                4: "ethnicity"})
        eigenstrat_result["ethnicity"] = eigenstrat_result[["ethnicity"]].replace([1,2,3], ["EU", "AF", "HISP"])
        return eigenstrat_result

    def PCAplot(csv_path, x, y, title,
                output_path = None, display = None, **kwargs):
        """normal PCA plot
        Args:
            csv_path: str, the EIGENSOFT smartPCA output file
            x, y: str, the variable(number of columns) to be used on x and y axis
            output_path: str, the saving folder(not the path of file) of the figure,
                         if provided, it will cover the previously provided path
                         JUST IN THIS RUN
            title: str, the title of the figure, and also the saving file name
                   of the figure
            display: Optional[bool], if provided, it will cover the previously
                    provided display setting JUST IN THIS RUN
            kwargs: any argument to be passed to sns.scatterplot method,
                    for plotting convenience, these setting will PERMANANTLY cover
                    the previous setting
        """
        data = self._read(path = csv_path)

        self._sns_setting = kwargs
        self._plot(data = data, x = x, y = y, title = title,
                   output_path = output_path, display = display)

    def PCAplot_one_vs_others(csv_path, x, y, one, title,
                              output_path = None, display = None, **kwargs):
        """split ethnicity into two parts, an individual ethnicity vs others

        Args:
            csv_path: str, the EIGENSOFT smartPCA output file
            x, y: str, the variable(number of columns) to be used on x and y axis
            one: the ethnicity to be chosen, all the other ethnicities specified
                 when Initialize will be consider as "others"(in figure legend
                 label "B & C & D" format)
            output_path: str, the saving folder(not the path of file) of the figure,
                         if provided, it will cover the previously provided path
                         JUST IN THIS RUN
            title: str, the title of the figure, and also the saving file name
                   of the figure
            display: Optional[bool], if provided, it will cover the previously
                    provided display setting JUST IN THIS RUN
            kwargs: any argument to be passed to sns.scatterplot method,
                    for plotting convenience, these setting will PERMANANTLY cover
                    the previous setting
        """
        data = self._read(path = csv_path)
        # create a copy of _ethnicities list, then generate list of "others" ethnicities
        others = self._ethnicities
        others.remove(one)
        # relabel the ethnicity column
        data["ethnicity"] = np.where(data["ethnicity"] == one.upper(),
                                     one.upper(),
                                     "&".join([ethnicity.upper() for ethnicity in others]))
        self._sns_setting = kwargs
        self._plot(data = data, x = x, y = y, title = title,
                   output_path = output_path, display = display)
