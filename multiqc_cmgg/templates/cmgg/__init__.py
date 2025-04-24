#!/usr/bin/env python
"""
=========
 CMGG
=========

CMGG Custom report template

"""
import os

template_dir = os.path.dirname(__file__)
template_parent = "default"
base_fn = "base.html"

import logging
logging.info("Custom template being loaded")
logging.info(f"path template{template_dir}")