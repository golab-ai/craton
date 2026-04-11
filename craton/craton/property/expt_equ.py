import os,sys
import math
#from ..database.mongodb import MongoDB
import json
import csv
from copy import deepcopy
import random
from pathlib import Path
from ..utils.figure import DrawFigure
from ..utils.numerical_algorithm import linear_fitting, rmse_calculate, rmse_r_calculate
import numpy as np


property_dict = {
	"density_of_liquid": "density_of_liquid",
	"density": "density_of_liquid",
	"density of liquid": "density_of_liquid",
	"liquid density": "density_of_liquid",
	"rho": "density_of_liquid",
	"vapor_pressure": "vapor_pressure",
	"vapor pressure": "vapor_pressure",
	"VP": "vapor_pressure",
	"coefficient_of_thermal_expansion_of_liquid": "coefficient_of_thermal_expansion_of_liquid",
	"Bliq": "coefficient_of_thermal_expansion_of_liquid",
	"diffusion_coefficient_in_air": "diffusion_coefficient_in_air",
	"diffusion coefficient in air": "diffusion_coefficient_in_air",
	"diffusion in air": "diffusion_coefficient_in_air",
	"air diffusion coefficient": "diffusion_coefficient_in_air",
	"air diffusion": "diffusion_coefficient_in_air",
	"diffusion_coefficient_at_infinite_dilution_in_water": "diffusion_coefficient_at_infinite_dilution_in_water",
	"diffusion coefficient at infinite dilution in water": "diffusion_coefficient_at_infinite_dilution_in_water",
	"diffusion_coefficient": "diffusion_coefficient_at_infinite_dilution_in_water",
	"diffusion coefficient in water": "diffusion_coefficient_at_infinite_dilution_in_water",
	"dilution diffusion coefficient": "diffusion_coefficient_at_infinite_dilution_in_water",
	"dilution diffusion coefficient in water": "diffusion_coefficient_at_infinite_dilution_in_water",
	"D0": "diffusion_coefficient_at_infinite_dilution_in_water",
	"explosion_limit_in_air": "explosion_limit_in_air",
	"explosion limit in air": "explosion_limit_in_air",
	"explosion": "explosion_limit_in_air",
	"explosion in air": "explosion_limit_in_air",
	"flash point": "explosion_limit_in_air",
	"autoigintion point": "explosion_limit_in_air",
	"autoigintion temperature": "explosion_limit_in_air",
	"threshold_limit_value": "threshold_limit_value",
	"heat_capacity_of_gas": "heat_capacity_of_gas",
	"heat capacity of gas": "heat_capacity_of_gas",
	"capacity of gas": "heat_capacity_of_gas",
	"gas heat capacity": "heat_capacity_of_gas",
	"gas capacity": "heat_capacity_of_gas",
	"heat_capacity_of_liquid": "heat_capacity_of_liquid",
	"heat_capacity_of_solid": "heat_capacity_of_solid",
	"heat capacity of solid": "heat_capacity_of_solid",
	"solid heat capacity": "heat_capacity_of_solid",
	"enthalpy_of_combustion": "enthalpy_of_combustion",
	"head of combustion": "enthalpy_of_combustion",
	"enthalpy of combustion": "enthalpy_of_combustion",
	"combustion heat": "enthalpy_of_combustion",
	"combustion enthalpy": "enthalpy_of_combustion",
	"enthalpy_of_fusion": "enthalpy_of_fusion",
	"heat of fusion": "enthalpy_of_fusion",
	"enthalpy of fusion": "enthalpy_of_fusion",
	"fusion heat": "enthalpy_of_fusion",
	"fusion enthalpy": "enthalpy_of_fusion",
	"enthalpy_of_vaporization": "enthalpy_of_vaporization",
	"enthalpy of vaporization": "enthalpy_of_vaporization",
	"heat of vaporization": "enthalpy_of_vaporization",
	"HOV": "enthalpy_of_vaporization",
	"vaporization enthalpy": "enthalpy_of_vaporization",
	"vaporization heat": "enthalpy_of_vaporization",
	"entropy_of_gas": "entropy_of_gas",
	"entropy of gas": "entropy_of_gas",
	"gas entropy": "entropy_of_gas",
	"entropy_of_formation": "entropy_of_formation",
	"entropy of formation": "entropy_of_formation",
	"formation entropy": "entropy_of_formation",
	"enthalpy_of_formation": "enthalpy_of_formation",
	"enthalpy of formation": "enthalpy_of_formation",
	"heat of formation": "enthalpy_of_formation",
	"formation enthalpy": "enthalpy_of_formation",
	"formation heat": "enthalpy_of_formation",
	"gibbs_energy_of_formation": "gibbs_energy_of_formation",
	"gibbs energy of formation": "gibbs_energy_of_formation",
	"gibbs free energy of formation": "gibbs_energy_of_formation",
	"free energy of formationformation free energy": "gibbs_energy_of_formation",
	"formation gibbs free energy": "gibbs_energy_of_formation",
	"formation gibbs energy": "gibbs_energy_of_formation",
	"helmholtz_energy_of_formation": "helmholtz_energy_of_formation",
	"helmholtz energy of formation": "helmholtz_energy_of_formation",
	"helmholtz free energy of formation": "helmholtz_energy_of_formation",
	"formation helmholtz energy": "helmholtz_energy_of_formation",
	"internal_energy_of_formation": "internal_energy_of_formation",
	"internal energy of formation": "internal_energy_of_formation",
	"energy of formation": "internal_energy_of_formation",
	"formation energy": "internal_energy_of_formation",
	"formation internal energy": "internal_energy_of_formation",
	"solubility_in_water": "solubility_in_water",
	"solubility": "solubility_in_water",
	"solubility in water": "solubility_in_water",
	"solubility_in_water_containing_salt": "solubility_in_water_containing_salt",
	"solubility in water containing salt": "solubility_in_water_containing_salt",
	"solubility in salt": "solubility_in_water_containing_salt",
	"solubility_of_gas_in_water": "solubility_of_gas_in_water",
	"solubility of gas in water": "solubility_of_gas_in_water",
	"solubility of gas": "solubility_of_gas_in_water",
	"gas solubility": "solubility_of_gas_in_water",
	"surface_tension": "surface_tension",
	"surface tension": "surface_tension",
	"thermal_conductivity_of_gas": "thermal_conductivity_of_gas",
	"thermal conductivity of gas": "thermal_conductivity_of_gas",
	"conductivity of gas": "thermal_conductivity_of_gas",
	"gas conductivity": "thermal_conductivity_of_gas",
	"gas thermal conductivity": "thermal_conductivity_of_gas",
	"thermal_conductivity_of_liquid": "thermal_conductivity_of_liquid",
	"thermal conductivity": "thermal_conductivity_of_liquid",
	"conductivity": "thermal_conductivity_of_liquid",
	"thermal conductivity of liquid": "thermal_conductivity_of_liquid",
	"conductivity of liquid": "thermal_conductivity_of_liquid",
	"liquid conductivity": "thermal_conductivity_of_liquid",
	"liquid thermal conductivity": "thermal_conductivity_of_liquid",
	"viscosity_of_gas": "viscosity_of_gas",
	"viscosity of gas": "viscosity_of_gas",
	"gas viscosity": "viscosity_of_gas",
	"viscosity_of_liquid": "viscosity_of_liquid",
	"viscosity": "viscosity_of_liquid",
	"viscosity of liquid": "viscosity_of_liquid",
	"liquid viscosity": "viscosity_of_liquid",
	"henrys_law_constant_for_compound_in_water": "henrys_law_constant_for_compound_in_water",
	"henrys law constant": "henrys_law_constant_for_compound_in_water",
	"henrys law constant for compound in water": "henrys_law_constant_for_compound_in_water",
	"henrys constant": "henrys_law_constant_for_compound_in_water",
	"henrys_law_constant_for_gas_in_water": "henrys_law_constant_for_gas_in_water",
	"henrys law constant for gas in water": "henrys_law_constant_for_gas_in_water",
	"henrys law constant for gas": "henrys_law_constant_for_gas_in_water",
	"gas henrys law constant": "henrys_law_constant_for_gas_in_water",
	"henrys constant for gas in water": "henrys_law_constant_for_gas_in_water",
	"henrys constant for gas": "henrys_law_constant_for_gas_in_water",
	"gas henrys constant": "henrys_law_constant_for_gas_in_water",
	"octanol_water_partition_coefficient": "octanol_water_partition_coefficient",
	"octanol water partition coefficient": "octanol_water_partition_coefficient",
	"partition coefficient": "octanol_water_partition_coefficient",
	"logP": "octanol_water_partition_coefficient",
	"logp": "octanol_water_partition_coefficient",
    "hydration_free_energy":"hydration_free_energy",
	"critical_temperature": "critical_temperature",
	"critical temperature": "critical_temperature",
	"crit. temp.": "critical_temperature",
	"Tc": "critical_temperature",
	"critical_pressure": "critical_pressure",
	"critical pressure": "critical_pressure",
	"crit. pressure": "critical_pressure",
	"Pc": "critical_pressure",
	"critical_volume": "critical_volume",
	"critical volume": "critical_volume",
	"crit. volume": "critical_volume",
	"Vc": "critical_volume",
	"critical_density": "critical_density",
	"critical density": "critical_density",
	"crit. density": "critical_density",
	"rhoc": "critical_density",
	"critical_compressibility": "critical_compressibility",
	"critical compressibility": "critical_compressibility",
	"crit. compressibility": "critical_compressibility",
	"Zc": "critical_compressibility",
	"acentric_factor": "acentric_factor",
	"acentric factor": "acentric_factor",
	"omega": "acentric_factor",
	"melting_point": "melting_point",
	"freezing point": "melting_point",
	"melting point": "melting_point",
	"boiling_point": "boiling_point",
	"boiling point": "boiling_point",
	"refractive_index": "refractive_index",
	"refractive index": "refractive_index",
	"solubility_parameter": "solubility_parameter",
	"solubility parameter": "solubility_parameter",
	"dipole_moment": "dipole_moment",
	"dipole moment": "dipole_moment",
	"dipole": "dipole_moment",
	"liquid_volume": "liquid_volume",
	"liquid volume": "liquid_volume",
	"van_der_waals_area": "van_der_waals_area",
	"vdw_area": "van_der_waals_area",
	"area": "van_der_waals_area",
	"vdw area": "van_der_waals_area",
	"van der waals area": "van_der_waals_area",
	"van_der_waals_volume": "van_der_waals_volume",
	"vdw_volume": "van_der_waals_volume",
	"vdw volume": "van_der_waals_volume",
	"van der waals volume": "van_der_waals_volume",
	"radius_of_gyration": "radius_of_gyration",
	"radius of gyration": "radius_of_gyration",
	"rg": "radius_of_gyration",
    'recommended_exposure_limit_mass_basis_NIOSH':'recommended_exposure_limit_mass_basis_NIOSH',
    'recommended_exposure_limit_vol_basis_NIOSH':'recommended_exposure_limit_vol_basis_NIOSH',
    'max_workplace_conc_vol_basis_germany':'max_workplace_conc_vol_basis_germany',
    'max_workplace_conc_mass_basis_germany':'max_workplace_conc_mass_basis_germany',
    'permissible_exposure_limit_mass_basis_OSHA':'permissible_exposure_limit_mass_basis_OSHA',
    'permissible_exposure_limit_vol_basis_OSHA':'permissible_exposure_limit_vol_basis_OSHA',
    'threshold_limit_vol_basis_ACGIH':'threshold_limit_vol_basis_ACGIH',
    'threshold_limit_mass_basis_ACGIH':'threshold_limit_mass_basis_ACGIH',
    'flash_point':'flash_point',
    'upper_explosive_limit':'upper_explosive_limit',
    'lower_explosive_limit':'lower_explosive_limit',
}

property_arr = [
    "density_of_liquid",    
    "melting_point",
    "boiling_point",
    "vapor_pressure",
    "VLE curve",
    "critical_temperature",
    "critical_pressure",
    "critical_density",
    "critical_volume",
    "critical_compressibility",
    "enthalpy_of_vaporization",
    "enthalpy_of_fusion",
    "enthalpy_of_combustion",
    "entropy_of_gas",
    "internal_energy_of_formation",
    "enthalpy_of_formation",
    "entropy_of_formation",
    "helmholtz_energy_of_formation",
    "gibbs_energy_of_formation",
    "solubility_of_gas_in_water",
    "solubility_in_water",
    "solubility_in_water_containing_salt",
    "solubility_parameter",
    "henrys_law_constant_for_gas_in_water",
    "henrys_law_constant_for_compound_in_water",
    "octanol_water_partition_coefficient",
    "hydration_free_energy",
    "pka",
    "heat_capacity_of_gas",
    "heat_capacity_of_liquid",
    "heat_capacity_of_solid",
    "coefficient_of_thermal_expansion_of_liquid",
    "isothermal compressibility",
    "isochoric heat capacity",
    "Joule-Thomson coefficient",
    "sound speed",
    "dielectron_constant (permittivity)",
    "surface_tension",
    "viscosity_of_liquid",
    "viscosity_of_gas",
    "thermal_conductivity_of_gas",
    "thermal_conductivity_of_liquid",
    "diffusion_coefficient_in_air",
    "diffusion_coefficient_at_infinite_dilution_in_water",
    "conductivity",
    "van_der_waals_volume",
    "van_der_waals_area",
    "radius_of_gyration",
    "dipole_moment",
    "acentric_factor",
    "refractive_index",
    "liquid_volume",
    "flash_point",
    "lower_explosive_limit",
    "upper_explosive_limit",
    "threshold_limit_vol_basis_ACGIH",
    "threshold_limit_mass_basis_ACGIH",
    "permissible_exposure_limit_vol_basis_OSHA",
    "permissible_exposure_limit_mass_basis_OSHA",
    "recommended_exposure_limit_vol_basis_NIOSH",
    "recommended_exposure_limit_mass_basis_NIOSH",
    "max_workplace_conc_vol_basis_germany",
    "max_workplace_conc_mass_basis_germany",
]

log_property = ["vapor_pressure",
                "solubility_of_gas_in_water",
                "solubility_in_water",
                "solubility_in_water_containing_salt",
                "henrys_law_constant_for_gas_in_water",
                "henrys_law_constant_for_compound_in_water",
                "van_der_waals_area",
            ]

def yaws_den_function(paras,T):
    A = paras[0]
    B = paras[1]
    n = paras[2]
    Tc = paras[3]
    return A*(B**(-1*(1 - T/Tc)**n))

def yaws_polynomial_function(paras,A):
    value = 0.0
    for ii,pp in enumerate(paras):
        value += pp*A**ii
    return value

def yaws_diffusion_coefficient_at_infinite_dilution_in_wate(paras,A):
    return 10**(paras[0] + paras[1]/A)

def yaws_entropy_function(paras,A):
    return paras[0] + paras[1] / A + paras[2] * A

def yaws_solubility_function(paras,A):
    return 10** (paras[0] + paras[1] / A + paras[2] / A / A)

def yaws_solubility_salt_function(paras,A):
    return 10 ** (paras[0] + paras[1] * A + paras[2] * A * A)

def yaws_exponential_function(paras,A):
    return paras[0] * (1 - A / paras[2]) ** paras[1]

def yaws_logarithmic_function(paras,A):
    value = paras[0] + paras[1] / A + paras[2]*math.log10(A)
    for ii,pp in enumerate(paras[3:]):
        value += pp * A ** (ii+1)
    value = 10**value
    return value

def yaws_viscosity_of_liquid_function(paras,A):
    return 10 ** (paras[0] + paras[1]/A + paras[2]*A + paras[3]*A*A)

def yaws_function(property,paras,RTmax,RTmin,Tmin,Tmax=None,Tval=None):
    __functions = {
                    "density_of_liquid":yaws_den_function,
                    "coefficient_of_thermal_expansion_of_liquid":yaws_exponential_function,
                    "diffusion_coefficient_in_air":yaws_polynomial_function,
                    "diffusion_coefficient_at_infinite_dilution_in_water":yaws_diffusion_coefficient_at_infinite_dilution_in_wate,
                    "enthalpy_of_vaporization":yaws_exponential_function,
                    "vapor_pressure":yaws_logarithmic_function,
                    "enthalpy_of_formation":yaws_polynomial_function,
                    "entropy_of_formation":yaws_entropy_function,
                    "entropy_of_gas":yaws_polynomial_function,
                    "gibbs_energy_of_formation":yaws_polynomial_function,
                    "heat_capacity_of_gas":yaws_polynomial_function,
                    "heat_capacity_of_liquid":yaws_polynomial_function,
                    "heat_capacity_of_solid":yaws_polynomial_function,
                    "helmholtz_energy_of_formation":yaws_polynomial_function,
                    "henrys_law_constant_for_gas_in_water":yaws_logarithmic_function,
                    "internal_energy_of_formation":yaws_polynomial_function,
                    "solubility_in_water":yaws_solubility_function,
                    "solubility_in_water_containing_salt":yaws_solubility_salt_function,
                    "solubility_of_gas_in_water":yaws_logarithmic_function,
                    "surface_tension":yaws_exponential_function,
                    "thermal_conductivity_of_gas":yaws_polynomial_function,
                    "thermal_conductivity_of_liquid":yaws_polynomial_function,
                    "viscosity_of_gas":yaws_polynomial_function,
                    "viscosity_of_liquid":yaws_viscosity_of_liquid_function,
                   }
    if property not in __functions:
        return "the property: %s is not be supported"%property
    property = property.lower()
    if Tmax is None:
        if RTmin > Tmin and RTmax < Tmin:
            return ["the temperature not in the range from %.2f to %.2f" %(RTmin,RTmax)]
        else:
            return __functions[property](paras,Tmin)
    else:
        if Tmin >= Tmax:
            return ["Tmin: %f lager than Tmax: %f"%(Tmin,Tmax)]
        if Tmin > RTmax or Tmax < RTmin:
            return ["the temperature range %f - %f not in the setting range %f - %f" %(Tmin,Tmax,RTmin,RTmax)]
        ss = ''
        if Tmin < RTmin:
            ss += "the input Tmin: %f lower than %.3f, so the Tmin is changed to %.3f\n"%(Tmin,RTmin,RTmin)
            Tmin = RTmin
        if Tmax > RTmax:
            ss += "the input Tmax: %f lager than %.3f, so the Tmax is changed to %.3f\n"%(Tmax,RTmax,RTmax)
            Tmax = RTmax
        values = [ss,]
        nn = int((Tmax - Tmin) / Tval)
        for ii in range(nn):
            T = Tmin + ii * Tval
            values.append([T, __functions[property](paras,T)])
        if Tmax == Tmin + nn*Tval:
            values.append([Tmax, __functions[property](paras,Tmax)])
        return values

def require_expt_data(molecules,properties,molecule_type=None,sources=None,temperatures=None,pressures=None,condinations=None,):
    from ..database.mongodb import MongoDB
    __function = {
        "yaws_function":yaws_function,
    }

    db = MongoDB()
    if molecule_type is None:
        molecule_type="inchi_key"
    if not isinstance(molecules,list):
        molecules = [molecules]
    if not isinstance(properties,list):
        properties = [properties]
    tmp_properties = {}
    for pp in properties:
        tmp_properties[property_dict[pp]] = pp
    properties = list(tmp_properties.keys())
    if sources is not None:
        if not isinstance(sources,list):
            sources = [sources]

    selector = {molecule_type:{"$in":molecules},"property":{"$in":properties}}
    if sources is not None:
        selector["source"]={"$in":sources}
    datas = []
    for doc in db.exptdata_coll.find(selector):
        datas.append(doc)
    results = []
    for doc in datas:
        if doc["assign"] == "function":
            if temperatures is not None and isinstance(temperatures,list):

                doc["values"] = __function[doc["function_method"]](doc["property"],
                                                                   doc["paras"],
                                                                   doc["Tmax"],
                                                                   doc["Tmin"],
                                                                   temperatures[0],
                                                                   Tmax=temperatures[1],
                                                                   Tval=temperatures[2]
                                                                   )
            else:
                T = temperatures if temperatures is not None else 298.15
                if T < doc["Tmin"]:
                    T = doc["Tmin"]
                elif T > doc["Tmax"]:
                    T = doc["Tmax"]
                #else:
                #    T = 298.15
                doc["values"] = ["",[T,__function[doc["function_method"]](doc["property"],doc["paras"],doc["Tmax"],doc["Tmin"],T)]]
            if doc["values"][0] == "":
                results.append({"error":None,"molecule_type":molecule_type,"molecule":doc[molecule_type],
                                "property":tmp_properties[doc['property']],"values":doc["values"][1:],"unit":doc["unit"]})

            else:
                results.append({"error":f"{molecule_type} {doc[molecule_type]} {tmp_properties[doc['property']]}: {doc['values'][0]}"})
        else:
            results.append({"error":None,"molecule_type":molecule_type,"molecule":doc[molecule_type],
                            "property":tmp_properties[doc['property']],"values":doc["value"],"unit":doc["unit"]})
    return results

def get_all_csv_property(inf,T,fn=None):
    datas = json.loads(open(inf).read())
    
    __function = {
        "yaws_function":yaws_function,
    }
    results = []
    error = []
    
    for doc in datas:
        CAS_number = doc["CAS_number"] if "CAS_number" in doc else "-"
        IUPAC_Name = doc["IUPAC_Name"] if "IUPAC_Name" in doc else "-"
        if doc["assign"] == "function":
            if "smiles" in doc:
                if doc["property"] == "solubility_in_water":
                    results.append([doc["property"],CAS_number,IUPAC_Name,doc["formula"],doc["smiles"],T,doc["value_25oC"]])
                elif doc["property"] == "solubility_in_water_containing_salt":
                    results.append([doc["property"],CAS_number,IUPAC_Name,doc["formula"],doc["smiles"],T,doc["value_25oC_34472ppm"]])
                else:
                    try:
                        if T >= float(doc["Tmin"]) and T <= float(doc["Tmax"]):
                            doc["value"] = __function[doc["function_method"]](doc["property"],
                            	                                       doc["paras"],
                                	                                   doc["Tmax"],
                                    	                               doc["Tmin"],
                                        	                           T,
                                            	                       )
                            results.append([doc["property"],CAS_number,IUPAC_Name,doc["formula"],doc["smiles"],T,doc["value"],doc["Tmax"],doc["Tmin"]])
                    except:
                        error.append(doc)
        else:
            if "smiles" in doc:
                if doc["property"] == "explosion_limit_in_air":
                    if doc["flash_point"]:
                        results.append(["flash_point",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["flash_point"]])
                    if doc["lower_explosive_limit"]:
                        results.append(["lower_explosive_limit",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["lower_explosive_limit"]])
                    if doc["upper_explosive_limit"]:
                        results.append(["upper_explosive_limit",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["upper_explosive_limit"]])
                elif doc["property"] == "threshold_limit_value":
                    if doc["threshold_limit_vol_basis_ACGIH"]:
                        results.append(["threshold_limit_vol_basis_ACGIH",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["threshold_limit_vol_basis_ACGIH"]])
                    if doc["threshold_limit_mass_basis_ACGIH"]:
                        results.append(["threshold_limit_mass_basis_ACGIH",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["threshold_limit_mass_basis_ACGIH"]])
                    if doc["permissible_exposure_limit_vol_basis_OSHA"]:
                        results.append(["permissible_exposure_limit_vol_basis_OSHA",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["permissible_exposure_limit_vol_basis_OSHA"]])
                    if doc["permissible_exposure_limit_mass_basis_OSHA"]:
                        results.append(["permissible_exposure_limit_mass_basis_OSHA",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["permissible_exposure_limit_mass_basis_OSHA"]])
                    if doc["recommended_exposure_limit_vol_basis_NIOSH"]:
                        results.append(["recommended_exposure_limit_vol_basis_NIOSH",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["recommended_exposure_limit_vol_basis_NIOSH"]])
                    if doc["recommended_exposure_limit_mass_basis_NIOSH"]:
                        results.append(["recommended_exposure_limit_mass_basis_NIOSH",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["recommended_exposure_limit_mass_basis_NIOSH"]])
                    if doc["max_workplace_conc_vol_basis_germany"]:
                        results.append(["max_workplace_conc_vol_basis_germany",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["max_workplace_conc_vol_basis_germany"]])
                    if doc["max_workplace_conc_mass_basis_germany"]:
                        results.append(["max_workplace_conc_mass_basis_germany",doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["max_workplace_conc_mass_basis_germany"]])
                else:
                    if doc["value"]:
                        results.append([doc["property"],doc["CAS_number"],doc["IUPAC_Name"],doc["formula"],doc["smiles"],"-",doc["value"]])
                        
    if fn is not None:
        with open(fn,'w') as outf:
            writer = csv.writer(outf)
            writer.writerows(results)
        properties = set([rr[0] for rr in results])
        for prop in properties:
                this_data = [rr for rr in results if rr[0] == prop]		
                print(prop,len(this_data))
                       	
        
        
    return results,error

def get_isothermo_property(inf,fn=None):
    datas  = json.loads(open(inf).read())
    __function = {
        "yaws_function":yaws_function,
    }
    results = []
    error = []

    def get_T_range(Tmin,Tmax):
        for kk in range(int(Tmin),int(Tmin)+20,1):
            if kk % 10 == 0:
                break
        for mm in range(int(Tmax),int(Tmax)-20,-1):
            if mm % 10 == 0:
                break
        return kk,mm

    for doc in datas:
        if doc["assign"] == "function":
            CAS_number = doc["CAS_number"] if "CAS_number" in doc else "-"
            IUPAC_Name = doc["IUPAC_Name"] if "IUPAC_Name" in doc else "-"
            if "smiles" in doc:
                try:
                    Tmin = float(doc["Tmin"])
                    Tmax = float(doc["Tmax"])
                    Rmin,Rmax = get_T_range(Tmin,Tmax)
                    for T in range(Rmin,Rmax+5,10):
                        try:
                            _doc = deepcopy(doc)
                            _doc["value"] = __function[_doc["function_method"]](_doc["property"],
                        	                                       _doc["paras"],
                            	                                   _doc["Tmax"],
                                	                               _doc["Tmin"],
                                    	                           T,
                                        	                       )
                            _doc["temperature"] = T
                            results.append(_doc)
                        except:
                            error.append(doc)
                except:
                    error.append(doc)

    if fn is not None:
        with open(fn,'w') as outf:
            outf.write(json.dumps(results))
        properties = set([rr["property"] for rr in results])
        for prop in properties:
                this_data = [rr for rr in results if rr["property"] == prop]		
                print(prop,len(this_data))
                       	
        
        
    return results,error

def get_all_json_property(inf,T,fn=None):
    datas  = json.loads(open(inf).read())
    __function = {
        "yaws_function":yaws_function,
    }
    results = []
    error = []
    
    __function = {
        "yaws_function":yaws_function,
    }
    results = []
    error = []
    def _get_property(doc,pp):
        this_doc = deepcopy(doc)
        this_doc["property"] = pp
        this_doc["temperature"] = "-"
        this_doc["value"] = doc[pp]
        this_doc["Tmin"] = "-"
        this_doc["Tmax"] = "-"
        return this_doc
    
    
    for doc in datas:
        CAS_number = doc["CAS_number"] if "CAS_number" in doc else "-"
        IUPAC_Name = doc["IUPAC_Name"] if "IUPAC_Name" in doc else "-"
        if doc["assign"] == "function":
            if "smiles" in doc:
                doc["temperature"] = T
                if doc["property"] == "solubility_in_water":
                    doc["value"] = doc["value_25oC"]
                    doc["Tmin"] = "-"
                    doc["Tmax"] = "-"
                    results.append(doc)
                elif doc["property"] == "solubility_in_water_containing_salt":
                    doc["value"] = doc["value_25oC_34472ppm"]
                    doc["Tmin"] = "-"
                    doc["Tmax"] = "-"
                    results.append(doc)
                else:
                    try:
                        if T >= float(doc["Tmin"]) and T <= float(doc["Tmax"]):
                            doc["value"] = __function[doc["function_method"]](doc["property"],
                            	                                       doc["paras"],
                                	                                   doc["Tmax"],
                                    	                               doc["Tmin"],
                                        	                           T,
                                            	                       )
                            results.append(doc)
                    except:
                        error.append(doc)
        else:
            if "smiles" in doc:
                if doc["property"] == "explosion_limit_in_air":
                    if doc["flash_point"]:
                        results.append(_get_property(doc,"flash_point"))
                    if doc["lower_explosive_limit"]:
                        results.append(_get_property(doc,"lower_explosive_limit"))
                    if doc["upper_explosive_limit"]:
                        results.append(_get_property(doc,"upper_explosive_limit"))
                elif doc["property"] == "threshold_limit_value":
                    if doc["threshold_limit_vol_basis_ACGIH"]:
                        results.append(_get_property(doc,"threshold_limit_vol_basis_ACGIH"))
                    if doc["threshold_limit_mass_basis_ACGIH"]:
                        results.append(_get_property(doc,"threshold_limit_mass_basis_ACGIH"))
                    if doc["permissible_exposure_limit_vol_basis_OSHA"]:
                        results.append(_get_property(doc,"permissible_exposure_limit_vol_basis_OSHA"))
                    if doc["permissible_exposure_limit_mass_basis_OSHA"]:
                        results.append(_get_property(doc,"permissible_exposure_limit_mass_basis_OSHA"))
                    if doc["recommended_exposure_limit_vol_basis_NIOSH"]:
                        results.append(_get_property(doc,"recommended_exposure_limit_vol_basis_NIOSH"))
                    if doc["recommended_exposure_limit_mass_basis_NIOSH"]:
                        results.append(_get_property(doc,"recommended_exposure_limit_mass_basis_NIOSH"))
                    if doc["max_workplace_conc_vol_basis_germany"]:
                        results.append(_get_property(doc,"max_workplace_conc_vol_basis_germany"))
                    if doc["max_workplace_conc_mass_basis_germany"]:
                        results.append(_get_property(doc,"max_workplace_conc_mass_basis_germany"))
                else:
                    if doc["value"]:
                        doc["temperature"] = "-"
                        doc["Tmin"] = "-"
                        doc["Tmax"] = "-"
                        results.append(doc)
                        
    if fn is not None:
        with open(fn,'w') as outf:
            outf.write(json.dumps(results))
        properties = set([rr["property"] for rr in results])
        for prop in properties:
                this_data = [rr for rr in results if rr["property"] == prop]		
                print(prop,len(this_data))
                       	
        
        
    return results,error

class SplitDataSet:
    def __init__(self,inf,prop,style,value=None,output_dir=".",test_flag=False) -> None:
        self.molecules = json.loads(open(inf).read())
        self.style = style
        self.prop = prop
        self.output_dir = output_dir
        Path(self.output_dir).mkdir(exist_ok=True)
        self.value = value
        self.test_flag = test_flag
        
        if self.prop == "all":
            self.prop = list(set([mol["property"] for mol in self.molecules]))
        if not isinstance(self.prop,list):
            self.prop = [self.prop]
        normal_prop = []
        for pp in self.prop:
            normal_prop.append(property_dict[pp])
        self.prop = normal_prop
        self.molecules_dict = {pp:[] for pp in self.prop}
        for mol in self.molecules:
            if mol["property"] in self.prop and "unavil" not in mol:
                if mol["property"] in log_property:
                    mol["value"] = math.log10(mol["value"])
                self.molecules_dict[mol["property"]].append(mol)
    
    
    def run(self):
        __FUN = {
                "random": self.random_split,
                "arom": self.arom_split,
                "ring": self.ring_split,
                "halogen": self.halogen_split,
                "ha": self.ha_split,
                "zelement":self.zelement_split,
                "element":self.element_split,
                "select": self.select_split,
                "total": self.total_split,
            }
        for pp in self.prop:
            test_set,train_set,fn = __FUN["total"](pp)
            self.write_file(test_set,train_set,fn)
            print(pp,"total",len(train_set),len(test_set))
            

            test_set,train_set,fn = __FUN["random"](pp)
            self.write_file(test_set,train_set,fn)
            print(pp,"random",len(train_set),len(test_set))
            
            if not self.test_flag:
                continue
            
            test_set,train_set,fn = __FUN["arom"](pp)
            self.write_file(test_set,train_set,fn)
            print(pp,"arom",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["arom"](pp,ratio=0.1,pre_fn="_10")
            self.write_file(test_set,train_set,fn)
            print(pp,"arom_10",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["halogen"](pp)
            self.write_file(test_set,train_set,fn)
            print(pp,"halogen",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["halogen"](pp,ratio=0.1,pre_fn="_10")
            self.write_file(test_set,train_set,fn)
            print(pp,"halogen_10",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["ha"](pp,pre_fn="_20")
            self.write_file(test_set,train_set,fn)
            print(pp,"ha_20",len(train_set),len(test_set))
            
            self.value = 15
            test_set,train_set,fn = __FUN["ha"](pp,pre_fn="_15")
            self.write_file(test_set,train_set,fn)
            print(pp,"ha_15",len(train_set),len(test_set))
            
            self.value = ["ZCON","ZCOS","ZCOP","ZCNS","ZCONS","ZCONP","ZCOSP","ZCONSP"]
            test_set,train_set,fn = __FUN["zelement"](pp,pre_fn="_main")
            self.write_file(test_set,train_set,fn)
            print(pp,"zelement_main",len(train_set),len(test_set))
            
            #self.value = ["ZCOP","ZCNS","ZCONS","ZCONP","ZCOSP","ZCONSP"]
            #test_set,train_set,fn = __FUN["zelement"](pp,pre_fn="_main2")
            #self.write_file(test_set,train_set,fn)
            #print(pp,"zelement_main2",len(train_set),len(test_set))
            
            self.value = ["XFCl","ZFBr","XClBr","XFClBr"]
            test_set,train_set,fn = __FUN["zelement"](pp,kk=1,pre_fn="_halogen")
            self.write_file(test_set,train_set,fn)
            print(pp,"zelement_halogen",len(train_set),len(test_set))
            
            self.value = ["S"]
            test_set,train_set,fn = __FUN["element"](pp,pre_fn="_S")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_S",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["element"](pp,ratio=0.1,pre_fn="_S10")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_S10",len(train_set),len(test_set))
            
            self.value = ["Cl"]
            test_set,train_set,fn = __FUN["element"](pp,pre_fn="_Cl")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_Cl",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["element"](pp,ratio=0.1,pre_fn="_Cl10")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_Cl10",len(train_set),len(test_set))
            
            self.value = ["F"]
            test_set,train_set,fn = __FUN["element"](pp,pre_fn="_F")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_F",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["element"](pp,ratio=0.1,pre_fn="_F10")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_F10",len(train_set),len(test_set))
            
            self.value = ["Br","I"]
            test_set,train_set,fn = __FUN["element"](pp,pre_fn="_BrI")
            self.write_file(test_set,train_set,fn)
            print(pp,"element_BrI",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["ring"](pp,pre_fn="_nonar")
            self.write_file(test_set,train_set,fn)
            print(pp,"ring_nonar",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["select"](pp,ratio=0.1,pre_fn="_10")
            self.write_file(test_set,train_set,fn)
            print(pp,"select_10",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["select"](pp,ratio=0.3,pre_fn="_30")
            self.write_file(test_set,train_set,fn)
            print(pp,"select_30",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["select"](pp,ratio=0.5,pre_fn="_50")
            self.write_file(test_set,train_set,fn)
            print(pp,"select_50",len(train_set),len(test_set))
            
            test_set,train_set,fn = __FUN["select"](pp,ratio=0.8,pre_fn="_80")
            self.write_file(test_set,train_set,fn)
            print(pp,"select_80",len(train_set),len(test_set))
            
    def write_file(self,test_set,train_set,fn):
        if len(test_set) > 0:
            with open(f"{self.output_dir}/{fn}_test_set.csv",'w') as outf:
                writer = csv.writer(outf)
                writer.writerows(test_set)
        with open(f"{self.output_dir}/{fn}_train_set.csv",'w') as outf:
            writer = csv.writer(outf)
            writer.writerows(train_set)
        
        if len(test_set) > 0:
            with open(f"{self.output_dir}/{fn}_test_set_ft.csv",'w') as outf:
                writer = csv.writer(outf)
                datas = [["SMILES","value"]] + [[rr[5],rr[6]] for rr in test_set]
                writer.writerows(datas)  
        with open(f"{self.output_dir}/{fn}_train_set_ft.csv",'w') as outf:
            writer = csv.writer(outf)
            datas = [["SMILES","value"]] + [[rr[5],rr[6]] for rr in train_set]
            writer.writerows(datas)

    def total_split(self,pp,pre_fn=""):
        molecules = self.molecules_dict[pp]
        train_set = []
        test_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            train_set.append(value)
        return test_set, train_set, f"{pp}_total{pre_fn}"

    def random_split(self,pp,pre_fn=""):
        molecules = self.molecules_dict[pp]
        nn = len(molecules)
        mm = int(nn*0.2)
        test_ndx = random.sample(range(0, nn), mm)
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            if ii in test_ndx:
                
                test_set.append(value)
            else:
                train_set.append(value)
        return test_set, train_set, f"{pp}_random{pre_fn}"
        
    def arom_split(self,pp,ratio=None,pre_fn=""):
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            if len(set(mol["ring_property"]).intersection({"ar1","ar2"})) > 0:
                test_set.append(value)
            else:
                train_set.append(value)
        if ratio is not None:
            nn = len(test_set)
            mm = int(nn*ratio)
            change_ndx = random.sample(range(0, nn), mm)
            _tmp_test_set = []
            for ii,value in enumerate(test_set):
                if ii in change_ndx:
                    train_set.append(value)
                else:
                    _tmp_test_set.append(value)
            test_set = _tmp_test_set
            
        return test_set,train_set,f"{pp}_arom{pre_fn}"
            
    def ring_split(self,pp,pre_fn=""):
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            tt = list(set(mol["ring_property"]))
            if len(tt) == 1 and tt[0] == "nonar":
                test_set.append(value)
            else:
                train_set.append(value)

        return test_set,train_set,f"{pp}_ring{pre_fn}"
            
    def halogen_split(self,pp,ratio=None,pre_fn=""):
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            zlabs = mol["zelement"].split("-")
            if len(zlabs[1]) == 1 and len(zlabs[2]) == 1:
                train_set.append(value)
            else:
                test_set.append(value)
        if ratio is not None:
            nn = len(test_set)
            mm = int(nn*ratio)
            change_ndx = random.sample(range(0, nn), mm)
            _tmp_test_set = []
            for ii,value in enumerate(test_set):
                if ii in change_ndx:
                    train_set.append(value)
                else:
                    _tmp_test_set.append(value)
            test_set = _tmp_test_set

        return test_set,train_set,f"{pp}_halogen{pre_fn}"
    
    def ha_split(self,pp,pre_fn=""):
        if self.value is None:
            vv = 20
        else:
            vv = int(self.value)
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            if mol["heavy_atoms"] > vv:
                test_set.append(value)
            else:
                train_set.append(value)

        return test_set,train_set,f"{pp}_ha{pre_fn}"

    def zelement_split(self,pp,kk=0,pre_fn=""):
        if self.value is None:
            targets = ["ZCONSP",]
        else:
            if isinstance(self.value,list):
                targets = self.value
            else:
                targets = self.value.split(":")
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            ss = mol["zelement"].split("-")
            if ss[kk] in targets:
                test_set.append(value)
            else:
                train_set.append(value)
        
        return test_set,train_set,f"{pp}_zelement{pre_fn}"   
    
    def element_split(self,pp,ratio=None,pre_fn=""):
        if self.value is None:
            targets = ["S"]
        else:
            if isinstance(self.value,list):
                targets = self.value
            else:
                targets = self.value.split(":")
        molecules = self.molecules_dict[pp]
        test_set = []
        train_set = []
        for ii,mol in enumerate(molecules):
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            if len(set(targets).intersection(set(mol["element_count"].keys()))) > 0:
                test_set.append(value)
            else:
                train_set.append(value)
                
        if ratio is not None:
            nn = len(test_set)
            mm = int(nn*ratio)
            change_ndx = random.sample(range(0, nn), mm)
            _tmp_test_set = []
            for ii,value in enumerate(test_set):
                if ii in change_ndx:
                    train_set.append(value)
                else:
                    _tmp_test_set.append(value)
            test_set = _tmp_test_set            
    

        return test_set,train_set, f"{pp}_element{pre_fn}"
    
    def select_split(self,pp,ratio=0.1,pre_fn=""):
        molecules = self.molecules_dict[pp]
        
        elem_dict = {"arom":[],"P":[],"I":[],"Br":[],"S":[],"F":[],"Cl":[],"N":[],"O":[],"C":[],"H":[]}
        _order = ["P","I","Br","S","F","Cl","N","O","C"]
        for mol in molecules:
            value = [mol["property"],mol["inchi_key"],mol["IUPAC_Name"],mol["CAS_number"],mol["formula"],mol["smiles"],mol["value"],mol["temperature"],mol["Tmin"],mol["Tmax"]]
            if len(set(mol["ring_property"]).intersection({"ar1","ar2"})) > 0:
                elem_dict["arom"].append(value)
            else:
                for ee in _order:
                    if ee in mol["element_count"]:
                        elem_dict[ee].append(value)
                        break

        elem_dict["C"] = elem_dict["C"] + elem_dict["H"]
        del elem_dict["H"]
        train_set = []
        test_set = []
        for aa,bb in elem_dict.items():
            nn = len(bb)
            if nn > 0:
                if nn == 1:
                    train_set.append(bb[0])
                else:
                    mm = int(ratio * nn)
                    if mm == 0:
                        train_set.append(bb[0])
                        test_set.extend(bb[1:])
                    else:
                        change_ndx = random.sample(range(0, nn), mm)
                        for ii,mol in enumerate(bb):
                            if ii in change_ndx:
                                train_set.append(mol)
                            else:
                                test_set.append(mol)
        return test_set,train_set,f"{pp}_select{pre_fn}"
        

def get_figure(inpath):
    out_path = f"{inpath}/results"
    fig_path = f"{inpath}/results/figure"
    ffs = [ff for ff in os.listdir(out_path) if ff.find(".csv") != -1]
    for ff in ffs:
        print(ff)
        datas = csv.reader(open(f"{out_path}/{ff}"))
        arr = []
        for data in datas:
            arr.append(data)
        x = [float(rr[6]) for rr in arr[:-1] if rr[0][0] != "#"]
        y = [float(rr[7]) for rr in arr[:-1] if rr[0][0] != "#"]
        DrawFigure.diagonal_draw([x,y],name=ff[:-4],rmse=True,fitting_curve=True,save_path=fig_path,rrmse=True,data_nn=True)

def get_results(inpath,color_shift=0):
    expt_path = f"{inpath}/origin"
    pred_path = f"{inpath}/predict"
    out_path = f"{inpath}/results"
    fig_path = f"{inpath}/results/figure"
    
    datas = {vv:{} for vv in property_arr}
    ffs = [ff for ff in os.listdir(pred_path)]
    out_of_range = []
    
    for ff in ffs:
        print(ff)
        pref = "_".join(ff.split("_")[:-3])
        data_pred = csv.reader(open(f"{pred_path}/{ff}"))
        dict_pred = {rr[1]:rr[2] for rr in data_pred}
        
        data_expt = csv.reader(open(f"{expt_path}/{ff[:-7]}.csv"))
        dict_expt = {rr[5]:rr for rr in data_expt}
        
        compare = []
        for m,pv in dict_pred.items():
            if m in dict_expt:
                er = dict_expt[m]
                if er[2] not in ["1,cis-3-dimethyl-cis-2-ethylcyclopentane"]:
                    pv = float(pv)
                    ev = float(er[6])

                    ue = abs(ev-pv)
                    rue = abs((ev-pv)/ev)
                    compare.append([er[0],er[1],er[2],er[3],er[4],er[5],ev,pv,ue,rue,er[7],er[8],er[9]])
                    
        x = [rr[6] for rr in compare]
        y = [rr[7] for rr in compare]
        compare = sorted(compare,key=lambda x:x[8],reverse=True)
        out_of_range.extend(compare[:10])
        a,b = rmse_calculate(x,y)
        ra,rb = rmse_r_calculate(x,y)
        c,d,e = linear_fitting(x,y)
        compare.append([b,a,rb,ra,e,c,d])
        with open(f"{out_path}/{ff}",'w') as outf:
            writer = csv.writer(outf)
            writer.writerows(compare)
        DrawFigure.diagonal_draw([x,y],name=ff[:-4],rmse=True,fitting_curve=True,save_path=fig_path,rrmse=True,data_nn=False,color_shift=color_shift)
        this_pp = compare[0][0]
        if pref not in datas[this_pp]:
            datas[this_pp][pref] = {"train":None,"test":None}
        if ff.find("train") != -1:
            datas[this_pp][pref]["train"] = [b,a,rb,ra,e,c,d]
        else:
            datas[this_pp][pref]["test"] = [b,a,rb,ra,e,c,d]
        
    with open(f"{out_path}/out_of_range.csv",'w') as outf:
        writer = csv.writer(outf)
        writer.writerows(out_of_range)

    with open(f"{out_path}/total_results.txt",'w') as outf:
        for this_pp in property_arr:
            if len(datas[this_pp]) > 0:
                for aa,bb in datas[this_pp].items():
                    s = bb["test"]
                    t = bb["train"]
                    if t is None:
                        outf.write("%s train - - - - - - - test %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" %(aa,s[0],s[1],s[2],s[3],s[4],s[5],s[6]))
                    else:
                        outf.write("%s train %.4f %.4f %.4f %.4f %.4f %.4f %.4f" %(aa,t[0],t[1],t[2],t[3],t[4],t[5],t[6]))
                        if s is not None:
                            s = bb["test"]
                            outf.write(" test %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" %(s[0],s[1],s[2],s[3],s[4],s[5],s[6]))
                        else:
                            outf.write(" test - - - - - - -\n")
                            
def get_admet_results(inpath,color_shift=0):
    expt_path = Path(inpath,"input_data")
    pred_path = Path(inpath,"predict")
    out_path = Path(inpath,"results")
    fig_path = f"{inpath}/results/figure"
    
    ffs = [ff for ff in os.listdir(pred_path)]
    
    for ff in ffs:
        print(ff)
        data_pred = list(csv.reader(open(f"{pred_path}/{ff}")))
        data_expt = list(csv.reader(open(f"{expt_path}/{ff}")))
        dict_pred = {dd[1]:dd[2] for dd in data_pred[1:]}
        dict_expt = {dd[0]:dd[1] for dd in data_expt[1:]}
        
        compare = []
        for m,vv in dict_pred.items():
            if m in dict_expt:
                ev = float(dict_expt[m])
                pv = float(vv)
                
                ue = abs(ev-pv)
                if ev == 0.0:
                    rue = 0.1
                else:
                    rue = abs((ev-pv)/ev)
    
                compare.append([m,ev,pv,ue,rue])
        compare = sorted(compare,key=lambda x:x[3],reverse=True)
        x = [dd[1] for dd in compare]
        y = [dd[2] for dd in compare]
        print(sum([dd[3] for dd in compare])/len(x))
        with open(f"{out_path}/{ff}",'w') as outf:
            writer = csv.writer(outf)
            writer.writerows(compare)
        DrawFigure.diagonal_draw([x,y],name=ff[:-4],rmse=True,fitting_curve=True,save_path=fig_path,rrmse=True,data_nn=True,color_shift=color_shift)

def get_fine_tune_script(path):
    #ss1 = "python finetune.py --finetune --batch-size 64 --lr 3e-4 --input-file-name data/inputdata --output-file-name &TRAIN& \n"
    ss1 = "python finetune_unimol.py --finetune --batch-size 128 --lr 4e-4 --input-file-name data/inputdata --output-file-name &TRAIN&\n"
    #ss2 = "python finetune.py --pred-smiles --output-file-name &MODEL& --val-folder output_csv/ --val-name &PRED& \n"
    ss2 = "python finetune_unimol.py --pred-smiles --output-file-name &MODEL& --val-folder output_csv/ --val-name &PRED&\n"
    ffs= [ff[:-4] for ff in os.listdir(f"{path}/input_data") if ff.find("_train") != -1]
    sss = ""
    for ff in ffs:
        #test = ff[:-13] + "_test_set_ft"
        test = ff.replace("train","test")
        tt1 = ss1
        tt1 = tt1.replace("&TRAIN&",ff)
        sss += tt1
        tt1 = ss2
        tt1 = tt1.replace("&MODEL&",ff)
        tt1 = tt1.replace("&PRED&",ff)
        sss += tt1
        if Path(path,"input_data",f"{test}.csv").exists():
            tt1 = ss2
            tt1 = tt1.replace("&MODEL&",ff)
            tt1 = tt1.replace("&PRED&",test)
            sss += tt1
    with open(f"{path}/input_data/run.sh",'w') as outf:
        outf.write(sss)
    
def property_value_analyze(ff,output_dir=".",bin_num=10):
    def _get_bin_data(datas,bin_num):
        min_v = min([rr["value"] for rr in datas])
        max_v = max([rr["value"] for rr in datas])
        size = max_v - min_v
        val = size/bin_num
        labels = [min_v + val*ii for ii in range(bin_num)] + [max_v*1.05]
        
        x = [[] for __ in labels]
        for dd in datas:
            for ii in range(0,bin_num):
                if dd["value"] >= labels[ii] and dd["value"] < labels[ii+1]:
                    x[ii].append(dd)
                    break
        kk = 0
        mols = []
        for ii in range(bin_num,-1,-1):
            kk += len(x[ii])
            if kk > 20:
                break
            else:
                mols.extend(x[ii])
        
        kk = 0
        for ii in range(0,bin_num,1):
            kk += len(x[ii])
            if kk > 20:
                break
            else:
                mols.extend(x[ii])
                
        labels = [round(vv,3) for vv in labels]
        return labels[:-1],[len(xx) for xx in x[:-1]],mols
    
    out_of_range = {}
    js = json.loads(open(ff).read())
    prop_dict = {}
    for rr in js:
        if "unavil" not in rr:
            if rr["property"] not in prop_dict:
                prop_dict[rr["property"]] = []
            prop_dict[rr["property"]].append(rr)
    for aa,bb in prop_dict.items():
        labels,x,mols = _get_bin_data(bb,bin_num)
        out_of_range[aa] = mols
        DrawFigure.bar(x,name=aa,x_label=labels,save_path=output_dir)
    arrs = []
    keys = sorted(list(out_of_range.keys()))
    for aa in keys:
        bb = out_of_range[aa]
        for bbb in bb:
            arrs.append([bbb["property"],bbb["inchi_key"],bbb["IUPAC_Name"],bbb["CAS_number"],bbb["formula"],bbb["smiles"],bbb["value"]])
    with open(f"{output_dir}/out_of_range.csv",'w') as outf:
        writer = csv.writer(outf)
        writer.writerows(arrs)

def admet_value_analyze(ff,output_dir=".",bin_num=10):
    def _get_bin_data(datas,bin_num):
        min_v = min([float(rr[1]) for rr in datas])
        max_v = max([float(rr[1]) for rr in datas])
        size = max_v - min_v
        val = size/bin_num
        labels = [min_v + val*ii for ii in range(bin_num)] + [max_v*1.05]
        
        x = [[] for __ in labels]
        for dd in datas:
            for ii in range(0,bin_num):
                if float(dd[1]) >= labels[ii] and float(dd[1]) < labels[ii+1]:
                    x[ii].append(dd)
                    break
        kk = 0
        mols = []
        for ii in range(bin_num,-1,-1):
            kk += len(x[ii])
            if kk > 20:
                break
            else:
                mols.extend(x[ii])
        
        kk = 0
        for ii in range(0,bin_num,1):
            kk += len(x[ii])
            if kk > 20:
                break
            else:
                mols.extend(x[ii])
                
        labels = [round(vv,3) for vv in labels]
        return labels[:-1],[len(xx) for xx in x[:-1]],mols
    
    out_of_range = {}
    js = json.loads(open(ff).read())
    for aa,bb in js.items():
        labels,x,mols = _get_bin_data(bb,bin_num)
        out_of_range[aa] = mols
        DrawFigure.bar(x,name=aa,x_label=labels,save_path=output_dir)
    arrs = []
    keys = sorted(list(out_of_range.keys()))
    for aa in keys:
        bb = out_of_range[aa]
        for bbb in bb:
            arrs.append(bbb)
    with open(f"{output_dir}/out_of_range.csv",'w') as outf:
        writer = csv.writer(outf)
        writer.writerows(arrs)


def unavil_record_of_property(rf,pf,outf):
    js = json.loads(open(pf).read())
    dict_r = {key:[] for key in set([rr["property"] for rr in js])}
    unavils = csv.reader(open(rf))
    for rr in unavils:
        if rr[0][0] != "#":
            dict_r[rr[0]].append(rr[2])
    
    for rr in js:
        if rr["IUPAC_Name"] in dict_r[rr["property"]]:
            rr["unavil"] = 0 
    
    with open(outf,'w') as outf:
        outf.write(json.dumps(js))
    
def add_expt_data(inf,propf,outf):
    prop_datas = json.loads(open(propf).read())
    attrs = ['IUPAC_Name', 'CAS_number', 'formula', 'source', 'assign', 'property', 'value', 'unit','smiles']
    if inf[-4:] == ".csv":
        datas =[]
        _tmp = list(csv.reader(open(inf)))
        kks = _tmp[0]
        for rr in _tmp[1:]:
            datas.append({kk:rr[ii].strip() for ii,kk in enumerate(kks)})
        #for data in datas:
        #    data["value"] = float(data["value"])
        #    data["value_std"] = float(data["value_std"])
        #    data["calc_value_std"] = float(data["calc_value_std"])
        #    data["calc_value"] = float(data["calc_value"])
            
    elif inf[-4:] == "json":
        datas = json.loads(open(inf).read())
    
    for data in datas:
        for attr in attrs:
            if attr not in data:
                data[attr] = "-"
        data["chem_info_flag"] = False
    
    prop_datas.extend(datas)
    
    with open(outf,'w') as outf:
        outf.write(json.dumps(prop_datas))

def get_expt_info_file(chemf,propf,outf):
    js=json.loads(open(propf).read())
    js1 = json.loads(open(chemf).read()) 
    for rr in js:
        if "smiles" in rr:
            if "chem_info_flag" in rr and not rr["chem_info_flag"]:
                rr.update(js1[rr["IUPAC_Name"]])
                rr["smiles"] = js1[rr["IUPAC_Name"]]["smiles"]
                rr["inchi_key"] = js1[rr["IUPAC_Name"]]["inchi_key"]
                rr["chem_info_flag"] = True
    with open(outf,'w') as outf:
        outf.write(json.dumps(js))
        
def get_npy_data(ff):
    outf =  open(f"./total_results.txt",'w')
    data = np.load(ff,allow_pickle=True)
    arrs = data.tolist()
    for aa,bb in arrs.items():
        compare = {}
        pref = "_".join(ff.split("_")[:-3])
        name = aa.split("/")[-1][:-4]
        
        x = bb["train_tgt_arr"]
        y = bb["train_pred"]
    
        DrawFigure.diagonal_draw([x,y],name=name,rmse=True,fitting_curve=True,save_path=".",rrmse=True,data_nn=True,color_shift=1)
        a,b = rmse_calculate(x,y)
        ra,rb = rmse_r_calculate(x,y)
        c,d,e = linear_fitting(x,y)
        outf.write("%s train %.3f %.3f %.3f %.3f %.3f %.3f %.3f "%(ff,b,a,rb,ra,e,c,d))
        
        
        name = name.replace("train_set","test_set")
        x = bb["test_tgt_arr"]
        y = bb["test_pred"]
    
        DrawFigure.diagonal_draw([x,y],name=name,rmse=True,fitting_curve=True,save_path=".",rrmse=True,data_nn=True,color_shift=1)
        
        a,b = rmse_calculate(x,y)
        ra,rb = rmse_r_calculate(x,y)
        c,d,e = linear_fitting(x,y)
        outf.write("test %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n"%(b,a,rb,ra,e,c,d))
        

        

    