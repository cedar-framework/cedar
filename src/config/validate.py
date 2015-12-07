#!/usr/bin/env python2

import sys

from jsonschema import Draft4Validator, validate
import json


def include_fileref(parent, mykey,  dic):
    for key in dic.keys():
        if key == '$ref':
            fpath = '../' + dic[key][5:]
            schemafile = open(fpath, 'r')
            schemastr = schemafile.read()
            schemafile.close()
            config_schema = json.loads(schemastr)
            parent[mykey] = config_schema

        elif type(dic[key]) is dict:
            include_fileref(dic, key, dic[key])

def main():
    schemafile = open('schema.json', 'r')
    schemastr = schemafile.read()
    schemafile.close()
    config_schema = json.loads(schemastr)
    include_fileref(None, None, config_schema)
    Draft4Validator.check_schema(config_schema)
    configfile = open(sys.argv[1], 'r')
    configstr = configfile.read()
    configfile.close()
    config = json.loads(configstr)
    validate(config, config_schema)

if __name__ == '__main__':
    main()
