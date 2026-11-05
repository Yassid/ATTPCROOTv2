# Embeds the trigger report figures into report.html.in. Usage: build_report.py <output.html>
import base64, os, re, sys
H = os.path.dirname(os.path.abspath(__file__)); F = '/mnt/f/a1975/sim_rebuild/trig/report/'
s = open(H + '/report.html.in').read()
s = re.sub(r'\{\{(\w+)\}\}', lambda m: 'data:image/png;base64,' + base64.b64encode(open(F + m[1] + '.png', 'rb').read()).decode(), s)
open(sys.argv[1], 'w').write(s); print(len(s))
