import re
import numpy as np
from datetime import datetime, timedelta
from typing import List, Dict, Any
from txt_parsing import strip_invisibles
def parse_head_movements(file_path, patient_id):
    quaternion_list = []
    times_list = []
    

    