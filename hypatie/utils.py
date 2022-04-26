import os
from datetime import datetime

def is_recent(fname):
    lm = os.path.getmtime(fname)
    lm = datetime.utcfromtimestamp(lm)
    now = datetime.utcnow()
    dt = (now - lm).total_seconds()
    recent = True if dt < 86400 else False
    return recent
