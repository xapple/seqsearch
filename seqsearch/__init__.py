#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Special variables #
__version__ = '2.2.1'

# After sh==1.14.3 the object returned changed #
import sh
sh_version = int(sh.__version__.split('.')[0])
if sh_version > 1: sh = sh.bake(_return_cmd=True)