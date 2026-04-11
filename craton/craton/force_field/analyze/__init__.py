from copy import deepcopy
from .mmqm_analyze import qm_mm_analyze
from ...utils.figure import DrawFigure as DF
     
class FittingAnlayze:
    def __init__(self) -> None:
        pass
    @staticmethod
    def analyze_qm_mm(
            molecules,
            results_path="./",
            force_field=None,
            atom_type_file=None,
            optimizer="openmm",
            done_fitting=None,
            init_this_ff=None,
            validation_terms=None,
            optimize_flag=True,
        ):
        if validation_terms is None:
            validation_terms = ["energy", "pes", "esp_charge", "Bonds", "Angles", "Dihedrals", "rmsd", "Pair1n","hessian","freq"]
        return qm_mm_analyze(
            molecules,
            results_path,
            force_field=force_field,
            atom_type_file=atom_type_file,
            optimizer=optimizer,
            done_fitting=done_fitting,
            init_this_ff=init_this_ff,
            validation_terms=validation_terms,
            optimize_flag=optimize_flag
        )

    @staticmethod
    def show_figure_fitting(datas,save_path):
        pes_args = {"labels":["QM","MM"],"rmse":True,"fitting_curve":True,"xylabels":["Dihedral","Energy"],"save_path":save_path}  
        args = {"rmse":True,"fitting_curve":True,"xylabels":["QM","MM"],"save_path":save_path}
        #plt_figure = plot_figure.DrawFigure(save_path=save_path)
        figure_path = []
        for typ, data in datas.items():
            if typ in ["rmsd"]:
                continue
            if typ != "pes":
                if len(data["total"][0]) > 0:
                    this_args = deepcopy(args)
                    this_args["name"] = typ
                    figure_path.append(DF.diagonal_draw([[data["total"][0], data["total"][1]]], **this_args))
            else:
                for mol_name, d in data.items():
                    for dih, v in d.items():
                        ####unique_dihedral_para = v[3][1]
                        name = "%s_%s" % (mol_name, "-".join(map(lambda x: str(int(x) + 1), dih.split("-"))))
                        if len(v) >= 4:
                            if all(i == 0 for i in v[3][2]):
                                name = "_nonfit_" + name
                        this_args = deepcopy(pes_args)
                        this_args["name"] = name
                        figure_path.append(DF.pes1d_draw([v[0], v[1]], **this_args))
        return figure_path

    @staticmethod
    def show_figure_parameter(para_valid_data,save_path):
        span_dict = {
            # "Bonds": 0.08,  # Angstrom
            # "Angles": 24,  # degree
            # "Impropers": 80,  # degree
            "Bonds": 0.15,  # Angstrom
            "Angles": 40,  # degree
            "Impropers": 120,  # degree
        }
        #draw = plot_figure.DrawFigure(save_path=save_path)
        figure_path = []
        for typ, data in para_valid_data.items():
            if typ.startswith("Dihedral"):
                continue
            for term, d in data.items():
                if len(d[0]) > 1:
                    name = typ + "--" + term
                    d_list_qm, d_list_mm, is_fitting = d
                    figure_path.append(DF.violin([d_list_qm, d_list_mm], name=name, span=span_dict[typ],save_path=save_path))
        return figure_path
