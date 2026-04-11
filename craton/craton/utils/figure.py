import statistics

import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from .. import CRATON_CONFIGURE
tmp_path = CRATON_CONFIGURE["path"]["tmp"]

from .numerical_algorithm import rmse_calculate, linear_fitting,rmse_r_calculate

plt.rcParams.update({"font.size": 13})

class DrawFigure:
    def __init__(self, save_path="./"):
        self.save_path = save_path
    _color_ = ["k","b","r","g","y","c","m",
               "grey","navy","brown","sage","orange","teal","purple",
               "silver","skyblue","peru","lime","gold","cyan","pink"]
    _line_ = ["-","--","-.",":"]

    _marker_ = ["o","*","+","x","v","1","s","p","h","D","|",
                "^","2",">","3","<","4","H","d","_",]

    @staticmethod
    def line_draw(XX, name=None, rmse=False, fitting_curve=False, labels=None, xylabels=None, save_path=None):
        if not isinstance(XX[0][0],list):
            XX = [XX]

        if save_path is None:
            save_path = tmp_path
        
        Z = [[],[]]
        for X in XX:
            Z[0].extend(X[0])
            Z[1].extend(X[1])

        min_v_x = min(Z[0])
        max_v_x = max(Z[0])
        min_v_y = min(Z[1])
        max_v_y = max(Z[1])
        if abs(max_v_x) >= abs(min_v_x):
            min_v_x = min_v_x - max_v_x * 0.05
            max_v_x = max_v_x + max_v_x * 0.05
        else:
            min_v_x = min_v_x + min_v_x * 0.05
            max_v_x = max_v_x - min_v_x * 0.05

        if abs(max_v_y) >= abs(min_v_y):
            min_v_y = min_v_y - max_v_y * 0.05
            max_v_y = max_v_y + max_v_y * 0.05
        else:
            min_v_y = min_v_y + min_v_y * 0.05
            max_v_y = max_v_y - min_v_y * 0.05

        
        plt.figure()

        if xylabels is not None:
            plt.xlabel(xylabels[0])
            plt.ylabel(xylabels[1])
        if name is not None:
            plt.title(name)
        
        plt.axis([min_v_x, max_v_x, min_v_y, max_v_y])
        
        texts = []

        for ii,X in enumerate(XX):
            args = {"color":DrawFigure._color_[ii % len(DrawFigure._color_)],
                "marker": DrawFigure._marker_[ii % len(DrawFigure._marker_)],
                "linestyle":"-",
                "markersize":5,
                }
            if labels is not None:
                args["label"] = labels[ii]

            plt.plot(X[0], X[1], **args)
            
        if labels is not None:
            plt.legend(loc=1)

        plt.tight_layout()
        plt.savefig(f"{save_path}/_line_{name}.png")
        plt.close()
        return f"{save_path}/_line_{name}.png"


    @staticmethod
    def diagonal_draw(XX, name=None, rmse=False, fitting_curve=False, labels=None, xylabels=None, save_path=None,rrmse=False,data_nn=False,color_shift=0):
        if not isinstance(XX[0][0],list):
            XX = [XX]

        if save_path is None:
            save_path = tmp_path
        
        Z = []
        for X in XX:
            Z.extend(X[0])
            Z.extend(X[1])

        min_v = min(Z)
        max_v = max(Z)
        if abs(max_v) >= abs(min_v):
            min_v = min_v - max_v * 0.05
            max_v = max_v + max_v * 0.05
        else:
            min_v = min_v + min_v * 0.05
            max_v = max_v - min_v * 0.05

        dia_arr = [[min_v, max_v], [min_v, max_v]]

        plt.figure()

        if xylabels is not None:
            plt.xlabel(xylabels[0])
            plt.ylabel(xylabels[1])
        if name is not None:
            plt.title(name)
        
        plt.axis([min_v, max_v, min_v, max_v])
        plt.plot(dia_arr[0], dia_arr[1], "k-", linewidth=1)
        
        texts = []

        for ii,X in enumerate(XX):
            args = {"color":DrawFigure._color_[ii % len(DrawFigure._color_) + color_shift],
                "marker": DrawFigure._marker_[ii % len(DrawFigure._marker_)],
                "linestyle":"",
                "markersize":5,
                }
            if labels is not None:
                args["label"] = labels[ii]

            plt.plot(X[0], X[1], **args)
            
            text = ""
            if data_nn:
                text += "data: %d\n" %len(X[1])
            if rmse:
                a, b = rmse_calculate(X[0], X[1])
                #text += "MAE=%.3f\nRMSE=%.3f\n"%(b,a)
                text += "MAE=%.3f\n" %b
            if rrmse:
                ra,rb = rmse_r_calculate(X[0],X[1])
                #text += f"RMAE=%.3f%%\nRRMSE=%.3f%%\n"%(rb*100,ra*100)
                text += f"MAE%%=%.3f%%\n" %(rb*100)
            if fitting_curve:
                c, d, e = linear_fitting(X[0], X[1])
                fit_curve = [[min_v, max_v], [c * min_v + d, c * max_v + d]]
                plt.plot(fit_curve[0], fit_curve[1], color = args["color"], linestyle = "--", linewidth=1)
                
                if d >=0:
                    #text += "y=%.3f*x+%.3f\nR2=%.3f" %(c,d,e)
                    text += "R2=%.3f\n" %e
                else:
                    #text += "y=%.3f*x-%.3f\nR2=%.3f" %(c,abs(d),e)
                    text += "R2=%.3f\n" %e
            if text != "":
                texts.append(text)
        if labels is not None:
            plt.legend(loc=1)

        if len(texts) > 0:
            if len(texts) == 1:
                this_text = texts[0]
            else:
                this_text = f"{labels[0]}\n{texts[0]}" if labels is not None else texts[0]
                for ii,text in enumerate(texts[1:]):
                    this_text += f"{labels[ii+1]}\n{text}" if labels is not None else f"\n{text}"
            plt.text(
                    0.05,
                    0.95,
                    this_text,
                    va="top",
                    fontsize=11,
                    transform=plt.gca().transAxes,
                    bbox=dict(boxstyle="round", facecolor="lavender", alpha=0.5),
                    )
        plt.tight_layout()
        plt.savefig(f"{save_path}/_diagonal_{name}.png")
        plt.close()
        return f"{save_path}/_diagonal_{name}.png"

    @staticmethod
    def pes1d_draw(XX, 
                   name=None, 
                   rmse=False, 
                   fitting_curve=False,
                   labels = None,
                   xylabels=None,
                   save_path=None,
                   ):
        if not isinstance(XX[0][0],list):
            XX = [XX]

        if save_path is None:
            save_path = "./"
        
        Z0 = []
        Z1 = []
        tmp = deepcopy(XX)
        XX = []
        for X in tmp:
            X = np.transpose(X)
            X = sorted(X, key=lambda x: x[0])
            X = np.transpose(X)
            XX.append(X)
            Z0.extend(X[0])
            Z1.extend(X[1])
        min_v_x = min(Z0)
        max_v_x = max(Z0)
        if abs(max_v_x) >= abs(min_v_x):
            min_v_x = min_v_x - max_v_x * 0.05
            max_v_x = max_v_x + max_v_x * 0.05
        else:
            min_v_x = min_v_x + min_v_x * 0.05
            max_v_x = max_v_x - min_v_x * 0.05
        min_v_y = min(Z1)
        max_v_y = max(Z1)
        if abs(max_v_y) >= abs(min_v_y):
            min_v_y = min_v_y - max_v_y * 0.05
            max_v_y = max_v_y + max_v_y * 0.05
        else:
            min_v_y = min_v_y + min_v_y * 0.05
            max_v_y = max_v_y - min_v_y * 0.05
        
        plt.figure()
        if xylabels is not None:
            if len(xylabels) == 1:
                xylabels.append(xylabels[0])
            plt.xlabel(xylabels[0])
            plt.ylabel(xylabels[1])
        if name is not None:
            plt.title(name)
        
        plt.axis([min_v_x, max_v_x, min_v_y, max_v_y])

        texts = []
        for ii,X in enumerate(XX):
            args = {"color":DrawFigure._color_[ii % len(DrawFigure._color_)],
                    "linestyle": "-",
                    "marker": "o",
                    "linewidth":1,
                    "markersize":5,
                    }
            if labels is not None:
                args["label"] = labels[ii]
            plt.plot(X[0],X[1],**args)
            if ii > 0:
                text = ""
                if rmse:
                    a, b = rmse_calculate(XX[0][1], X[1])
                    text += "MAE=%.3f\nRMSE=%.3f\n" %(b,a)
                if fitting_curve:
                    c, d, e = linear_fitting(XX[0][1], X[1])
                    if d >=0:
                        text += "y=%.3f*x+%.3f\nR2=%.3f" %(c,d,e)
                    else:
                        text += "y=%.3f*x-%.3f\nR2=%.3f" %(c,abs(d),e)
                if text != "":
                    texts.append(text)

        if labels is not None:
            plt.legend(loc=1)

        if len(texts) > 0:
            if len(texts) == 1:
                this_text = texts[0]
            else:
                this_text = f"{labels[1]}\n{texts[0]}" if labels is not None else texts[0]
                for ii,text in enumerate(texts[1:]):
                    this_text += f"{labels[ii+2]}\n{text}" if labels is not None else f"\n{text}"
            plt.text(
                    0.05,
                    0.95,
                    this_text,
                    va="top",
                    fontsize=11,
                    transform=plt.gca().transAxes,
                    bbox=dict(boxstyle="round", facecolor="lavender", alpha=0.5),
                    )
            
        plt.tight_layout()
        plt.savefig(f"{save_path}/{name}.png")
        plt.close()
        return f"{save_path}/{name}.png"

    @staticmethod
    def violin(XX, name=None, span=None,save_path=None):
        if not isinstance(XX[0],list):
            XX = [XX]
        if save_path is None:
            save_path = tmp_path
        fig, ax = plt.subplots()
        parts = ax.violinplot(XX, showmeans=True)
        parts["bodies"][0].set_color("olive")
        parts["bodies"][1].set_color("lime")

        if span is not None:
            #mean1 = statistics.mean(list1)
            #mean2 = statistics.mean(list2)
            #y_mid = (mean1 + mean2) / 2
            y_mid = sum([statistics.mean(X) for X in XX]) / len(XX)
            ax.set_ylim(y_mid - span / 2, y_mid + span / 2)

        plt.title(name.replace("--", "  ").replace("$", " "))
        plt.tight_layout()
        plt.savefig(f"{save_path}/_violin_{name}.png")
        plt.close()
        return f"{save_path}/_violin_{name}.png"

    def pie(self, param, param1, param2):
        pass

    @staticmethod
    def bar(Y,name=None,x_label=None,labels=None,show_value=True,save_path=None):
        if save_path is None:
            save_path = tmp_path
        if not isinstance(Y[0],list):
            Y = [Y]
            #x_label=[x_label]
            if labels is not None:
                labels = [labels]

        figname = f"{save_path}/{name}.png"

        #if x_label is None:
        #    x_label = [ii+1 for ii in range(len(Y[0]))]
        
        X = [i + 1 for i in range(len(Y[0]))]

        size_x = 8 
        size_y = 8
        #plt.figure(figsize=(size_x, size_y))
        plt.figure()

        width = 0.2
        for ii in range(0, len(Y)):
            shift = (ii + 0.5) * width
            if labels is not None:
                plt.bar([x + shift for x in X], Y[ii], width, label=labels[ii])
            else:
                plt.bar([x + shift for x in X], Y[ii], width)
            if show_value:
                for jj, vv in enumerate(Y[ii]):
                    plt.text(X[jj] + shift, vv / 2.0, "%d" % vv, ha="center", fontsize=10)

        if x_label is not None:
            plt.xticks(X, labels=x_label, rotation=45,fontsize=8)
        #plt.xticks(X,labels=x_label[0])
        plt.title(name)
        plt.grid(True, linestyle=":", color="r", alpha=0.6)
        plt.legend()
        plt.savefig(figname)
        plt.close()
        return figname