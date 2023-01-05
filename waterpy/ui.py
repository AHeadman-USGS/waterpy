"""Objects and functions that aid in the transition from UI to main"""

import os


class Inputs:
    """Likely a temporary place holder while we work out how to pass
    between pyQT and here """
    input_directory = 'dir'
    output_directory = 'dir'
    characteristics = 'file'
    parameters = 'file'
    timeseries = 'file'
    twi = 'file'
    pet = 'hamon'
    snow = 'option'
    karst = 'option'
    rand = 'option'
    matrices = 'option'


def ini_write(file, Inputs):
    """Build the modelconfig.ini file from user inputs"""
    config_file = open(file, "w")
    lines = []
    # inputs are going to be brain hurty.
    lines.extend(("[Inputs] \n",
                  "input_dir = {} \n".format(Inputs.input_directory),
                  "characteristics_basin_file = ${Inputs:input_dir}" + "\\{} \n".format(Inputs.characteristics),
                  "parameters_land_type_file = ${Inputs:input_dir}" + "\\{} \n".format(Inputs.parameters),
                  "timeseries_file = ${Inputs:input_dir}" + "\\{} \n".format(Inputs.timeseries),
                  "twi_file = ${Inputs:input_dir}\\" + "{} \n".format(Inputs.twi),
                  "\n[Outputs] \n",
                  "output_dir = {} \n".format(Inputs.output_directory),
                  "output_filename = output.csv \n",
                  "output_filename_saturation_deficit_locals = output_saturation_deficit_locals.csv \n",
                  "output_filename_unsaturated_zone_storages = output_unsaturated_zone_storages.csv \n",
                  "output_filename_root_zone_storages = output_root_zone_storages.csv \n",
                  "output_filename_evaporations = output_evaporations.csv \n",
                  "output_report = report.html \n",
                  "\n[Options] \n",
                  "option_pet = {} \n".format(Inputs.pet),
                  "option_snowmelt = {}\n".format(Inputs.snow),
                  "option_karst = {}\n".format(Inputs.karst),
                  "option_randomize_daily_to_hourly = {}\n".format(Inputs.rand),
                  "option_write_output_matrices = {}\n".format(Inputs.matrices)))

    for line in lines:
        config_file.write(line)
    config_file.close()

    return os.path.abspath(file)
