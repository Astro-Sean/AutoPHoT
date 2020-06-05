#!/usr/bin/python

'''

Package used for load yaml files and
dumping/ updating yaml files to control AutoPhot Input

Input Parameters:
    - fname = name of yaml file
    - dname = location of fname files, if None location will be direcotyr of call_yaml.py

'''


class yaml_syntax(object):

    def __init__(self,filepath = None,dict_name = None,syntax = None):

         self.filepath  = filepath
         self.dict_name = dict_name
         self.syntax = syntax

    def load_vars(self):

        import yaml
        import os

        if self.syntax != None:
            file_path = os.path.join(self.syntax['wdir'], self.filepath )
        else:
            file_path = self.filepath


        with open(file_path, 'r') as stream:
            var = yaml.load(stream, Loader=yaml.FullLoader)

        if self.dict_name != None:
            data = var[self.dict_name]
        else:
            data = var

        return data

    def update_vars(self,key,new_val):



        import ruamel
        import ruamel.yaml


        with open(self.filepath,'r') as yamlfile:
            doc = ruamel.yaml.load(yamlfile, Loader=ruamel.yaml.RoundTripLoader,preserve_quotes=True)

            doc[self.dict_name][str(key)] =  new_val

        with open(self.filepath,'w') as yamlfile:
            ruamel.yaml.dump(doc,yamlfile, Dumper=ruamel.yaml.RoundTripDumper)

    def new_vars(self,key,val):

        import ruamel
        import ruamel.yaml
        doc = {self.dict_name:{key:val}}

        with open(self.filepath,'r') as yamlfile:


            cur_yaml = ruamel.yaml.load(yamlfile, Loader=ruamel.yaml.RoundTripLoader,preserve_quotes=True)
            cur_yaml.update(doc)

        with open(self.filepath,'w') as yamlfile:
            ruamel.yaml.dump(cur_yaml,yamlfile, Dumper=ruamel.yaml.RoundTripDumper)


    def create_yaml(fname,data):
        import yaml
        import os

        target_name = fname
        if '.yml' not in fname:
            fname+='.yml'

        data_new  = {os.path.basename(target_name.replace('.yml','')):data}
        with open(fname, 'w') as outfile:
            yaml.dump(data_new, outfile, default_flow_style=False)



