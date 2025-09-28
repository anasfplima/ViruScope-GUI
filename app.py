from pathlib import Path
from Bio import SeqIO
from shiny import App, ui, render, reactive, session
import tempfile
import string
import random
import os
import pandas as pd
import numpy as np
import shutil 
import atexit
import asyncio
import matplotlib.pyplot as plt
from Conservation_Scores_V2 import viruscope
from arolit import arolit
from safedir_generator import safe_delete_folder
import utilities_fasta as utilities
from matplotlib_venn import venn2


# Track all created temp dirs so we can clean them up when app exits
all_temp_dirs = []

def make_temp_dir(base_dir):
    tmp = tempfile.mkdtemp(dir=base_dir)
    all_temp_dirs.append(tmp)
    print(f"[TempDir] Created: {tmp}")
    return tmp

def cleanup_all_temp_dirs():
    for tmp in all_temp_dirs:
        if os.path.exists(tmp):
            shutil.rmtree(tmp, ignore_errors=True)
            print(f"[TempDir] Deleted: {tmp}")



# ------------------------------------------------

atexit.register(cleanup_all_temp_dirs)

project_dir = os.path.dirname(os.path.abspath(__file__))
www_dir = Path(__file__).parent / "www"

app_ui = ui.page_navbar(
    ui.head_content(
        ui.tags.link(rel="icon", href="favicon.ico"),
        ui.tags.link(rel='javascript', href='controls.js'),
        ui.tags.link(rel="stylesheet", href="styles.css"),
        ui.tags.link(
            rel="stylesheet",
            href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.7.2/css/all.min.css",
        )
    ),
    ui.nav_spacer(),

    # Home
    ui.nav_panel("Home",
            ui.div(
                ui.div(
                    ui.HTML(
                        """
                        <div class="home-left">
                            <h1 class= "home-title">Welcome to Viruscope</h1>
                            <p class= "home-subtitle">A suite of tools for viral primer retrieval, generation and scoring.</p>
                        </div>
                        <div class="home-right">
                            <?xml version="1.0" encoding="UTF-8"?><svg class="rotate-center" id="b" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 781.62 781.62"><g id="c"><g><polyline points="154.32 534.19 116.67 394.99 160.91 494.09 154.32 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="154.32 534.19 160.91 494.09 223.98 480.44 154.32 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="160.91 494.09 216.59 380.49 223.98 480.44 160.91 494.09" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="216.59 380.49 222.87 348.91 232.47 343.71 216.59 380.49" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="222.87 348.91 154.32 255.79 232.47 343.71 222.87 348.91" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="253.27 633.14 154.32 534.19 271.47 605.38 253.27 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="255.61 168.16 253.27 156.84 284.53 168.54 255.61 168.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="255.61 168.16 284.53 168.54 303.3 254.81 255.61 168.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="223.98 480.44 216.59 380.49 302.43 493.08 223.98 480.44" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="216.59 380.49 232.47 343.71 345.83 380.29 216.59 380.49" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="303.3 254.81 284.53 168.54 348.35 207.8 303.3 254.81" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="302.43 493.08 216.59 380.49 345.83 380.29 302.43 493.08" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="302.43 493.08 345.83 380.29 367.36 401.12 302.43 493.08" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="232.47 343.71 303.3 254.81 345.83 380.29 232.47 343.71" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="303.3 254.81 348.35 207.8 408.6 253.28 303.3 254.81" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="345.83 380.29 303.3 254.81 408.6 253.28 345.83 380.29" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="470.65 207.24 480.54 152.84 519.95 179.18 470.65 207.24" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="470.65 207.24 519.95 179.18 525.59 190.53 470.65 207.24" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="519.95 179.18 480.54 152.84 531.67 156.84 519.95 179.18" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="531.67 633.14 460.09 571.94 579.96 583.21 531.67 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="477.36 394.63 470.65 207.24 586.09 232.04 477.36 394.63" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="477.36 394.63 586.09 232.04 641.88 326.66 477.36 394.63" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="586.09 232.04 630.62 255.79 641.88 326.66 586.09 232.04" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="641.88 326.66 630.62 255.79 668.27 394.99 641.88 326.66" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="191.12 226.26 253.27 156.84 255.61 168.16 191.12 226.26" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="154.32 255.79 191.12 226.26 232.47 343.71 154.32 255.79" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="232.47 343.71 191.12 226.26 303.3 254.81 232.47 343.71" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="191.12 226.26 255.61 168.16 303.3 254.81 191.12 226.26" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="191.12 226.26 154.32 255.79 253.27 156.84 191.12 226.26" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="470.65 207.24 525.59 190.53 586.09 232.04 470.65 207.24" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="271.47 605.38 326.42 527.7 369.41 616.11 271.47 605.38" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="253.27 633.14 371.6 650.86 392.47 670.79 253.27 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="160.91 494.09 116.67 394.99 216.59 380.49 160.91 494.09" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="134.23 363.16 116.67 394.99 154.32 255.79 134.23 363.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="116.67 394.99 134.23 363.16 216.59 380.49 116.67 394.99" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="216.59 380.49 134.23 363.16 222.87 348.91 216.59 380.49" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="134.23 363.16 154.32 255.79 222.87 348.91 134.23 363.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="477.36 394.63 641.88 326.66 668.27 394.99 477.36 394.63" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="399.76 595.51 326.42 527.7 422.64 552.77 399.76 595.51" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="422.64 552.77 326.42 527.7 429.32 438.86 422.64 552.77" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="326.42 527.7 302.43 493.08 429.32 438.86 326.42 527.7" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="369.41 616.11 326.42 527.7 399.76 595.51 369.41 616.11" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="408.6 253.28 470.65 207.24 477.36 394.63 408.6 253.28" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="408.6 253.28 348.35 207.8 470.65 207.24 408.6 253.28" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="432.99 648.64 460.09 571.94 531.67 633.14 432.99 648.64" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="392.47 670.79 432.99 648.64 531.67 633.14 392.47 670.79" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="253.27 633.14 271.47 605.38 281.66 622.26 253.27 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="281.66 622.26 369.41 616.11 371.6 650.86 281.66 622.26" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="281.66 622.26 271.47 605.38 369.41 616.11 281.66 622.26" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="253.27 633.14 281.66 622.26 371.6 650.86 253.27 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="531.67 633.14 579.96 583.21 630.62 534.19 531.67 633.14" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="636.84 479.86 477.36 394.63 668.27 394.99 636.84 479.86" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="630.62 534.19 636.84 479.86 668.27 394.99 630.62 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="154.32 534.19 223.98 480.44 250.19 502 154.32 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="154.32 534.19 250.19 502 271.47 605.38 154.32 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="250.19 502 223.98 480.44 302.43 493.08 250.19 502" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="271.47 605.38 250.19 502 326.42 527.7 271.47 605.38" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="250.19 502 302.43 493.08 326.42 527.7 250.19 502" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="367.36 401.12 345.83 380.29 390.42 405.16 367.36 401.12" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="302.43 493.08 390.42 405.16 429.32 438.86 302.43 493.08" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="302.43 493.08 367.36 401.12 390.42 405.16 302.43 493.08" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="429.32 438.86 390.42 405.16 477.36 394.63 429.32 438.86" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="390.42 405.16 345.83 380.29 408.6 253.28 390.42 405.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="390.42 405.16 408.6 253.28 477.36 394.63 390.42 405.16" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="422.64 552.77 429.32 438.86 460.09 571.94 422.64 552.77" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="525.59 190.53 519.95 179.18 556.92 182.53 525.59 190.53" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="519.95 179.18 531.67 156.84 556.92 182.53 519.95 179.18" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="525.59 190.53 556.92 182.53 586.09 232.04 525.59 190.53" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="586.09 232.04 556.92 182.53 630.62 255.79 586.09 232.04" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="556.92 182.53 531.67 156.84 630.62 255.79 556.92 182.53" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="399.76 595.51 422.64 552.77 460.09 571.94 399.76 595.51" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="432.99 648.64 399.76 595.51 460.09 571.94 432.99 648.64" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="460.09 571.94 429.32 438.86 569.19 524.02 460.09 571.94" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="460.09 571.94 569.19 524.02 579.96 583.21 460.09 571.94" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="429.32 438.86 477.36 394.63 569.19 524.02 429.32 438.86" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="579.96 583.21 569.19 524.02 630.62 534.19 579.96 583.21" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="630.62 534.19 569.19 524.02 636.84 479.86 630.62 534.19" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="569.19 524.02 477.36 394.63 636.84 479.86 569.19 524.02" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="348.35 207.8 284.53 168.54 392.47 119.19 348.35 207.8" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="284.53 168.54 253.27 156.84 392.47 119.19 284.53 168.54" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="348.35 207.8 392.47 119.19 470.65 207.24 348.35 207.8" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="470.65 207.24 392.47 119.19 480.54 152.84 470.65 207.24" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="480.54 152.84 392.47 119.19 531.67 156.84 480.54 152.84" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="371.6 650.86 369.41 616.11 385.72 630.46 371.6 650.86" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="385.72 630.46 369.41 616.11 399.76 595.51 385.72 630.46" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="371.6 650.86 385.72 630.46 392.47 670.79 371.6 650.86" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="392.47 670.79 385.72 630.46 432.99 648.64 392.47 670.79" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><polyline points="385.72 630.46 399.76 595.51 432.99 648.64 385.72 630.46" fill="none" stroke=""none"" stroke-miterlimit="10" stroke-width="2"/><circle cx="154.32" cy="534.19" r="7.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="116.67" cy="394.99" r="2.52" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="8.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="1.31" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r=".6" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="1.51" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="2.42" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="7.25" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="2.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="222.87" cy="348.91" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="5.44" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="222.87" cy="348.91" r="6.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r=".81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="8.56" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="222.87" cy="348.91" r="4.03" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="5.94" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="9.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="3.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="3.12" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="156.84" r="4.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="5.84" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="4.63" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="6.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="4.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="9.06" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="1.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="5.34" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="3.42" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="7.05" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="3.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="2.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="1.51" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="1.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="3.93" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="7.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="6.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="367.36" cy="401.12" r="9.77" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r=".91" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="4.63" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="3.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="4.53" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="7.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="1.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="5.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="1.91" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="8.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="4.73" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="7.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="3.93" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="480.54" cy="152.84" r="6.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="9.57" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r=".5" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="2.82" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="7.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="3.02" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="156.84" r="8.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="480.54" cy="152.84" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="5.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="579.96" cy="583.21" r="4.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="6.24" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="3.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="5.64" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="4.63" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r=".6" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="1.31" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="3.02" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="641.88" cy="326.66" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="7.35" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="8.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="641.88" cy="326.66" r=".6" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="255.79" r="6.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="2.72" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="641.88" cy="326.66" r="8.76" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="668.27" cy="394.99" r="6.34" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="255.79" r="4.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="641.88" cy="326.66" r="9.47" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="2.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="156.84" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="4.03" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="4.03" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="4.13" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="7.05" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="2.82" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="2.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="232.47" cy="343.71" r="7.35" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="303.3" cy="254.81" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="255.61" cy="168.16" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="1.31" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r=".81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="156.84" r="8.76" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="5.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="191.12" cy="226.26" r="2.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="2.52" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="2.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="2.82" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="2.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="4.13" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="3.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="3.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="2.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="6.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="7.45" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="4.33" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="116.67" cy="394.99" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="160.91" cy="494.09" r="9.16" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="116.67" cy="394.99" r="8.66" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="2.92" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="116.67" cy="394.99" r="8.86" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r=".81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="116.67" cy="394.99" r="1.21" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="222.87" cy="348.91" r="8.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="216.59" cy="380.49" r="6.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="222.87" cy="348.91" r="4.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="255.79" r="2.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="134.23" cy="363.16" r="6.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="1.61" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="668.27" cy="394.99" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="641.88" cy="326.66" r="5.24" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="9.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="3.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="7.05" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="8.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="7.35" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="1.31" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="1.51" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="3.12" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="4.63" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="4.63" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="5.34" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="9.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="2.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="2.82" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="7.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="9.47" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="8.06" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="4.73" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="3.52" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="2.72" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="8.76" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="9.57" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="3.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="7.25" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="6.04" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="7.35" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="6.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="9.77" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="2.92" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="7.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="7.95" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="8.86" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="8.86" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="1.61" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="4.33" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="7.05" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="281.66" cy="622.26" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="633.14" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="6.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="4.03" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="579.96" cy="583.21" r="7.25" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="633.14" r="7.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="636.84" cy="479.86" r="8.76" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="668.27" cy="394.99" r="3.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="6.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="636.84" cy="479.86" r="7.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="6.24" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="668.27" cy="394.99" r="2.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="636.84" cy="479.86" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="9.97" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="5.64" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="3.73" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="4.33" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="7.45" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="5.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="154.32" cy="534.19" r="3.52" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="3.02" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="6.24" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="223.98" cy="480.44" r="6.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="2.92" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="8.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="7.15" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="3.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="271.47" cy="605.38" r="7.45" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="8.86" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="326.42" cy="527.7" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="250.19" cy="502" r="7.35" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="367.36" cy="401.12" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="9.77" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="3.02" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="367.36" cy="401.12" r=".5" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="2.42" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="2.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r=".7" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="5.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="367.36" cy="401.12" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="302.43" cy="493.08" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="8.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="8.56" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="7.05" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="8.76" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="345.83" cy="380.29" r="8.56" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="4.13" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="1.31" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="5.24" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="408.6" cy="253.28" r="7.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="390.42" cy="405.16" r="7.45" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="6.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="7.25" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="2.92" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r="3.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="7.45" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="1.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r="8.56" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="3.73" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="8.16" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="156.84" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="519.95" cy="179.18" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r="4.23" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="3.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="525.59" cy="190.53" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="255.79" r="3.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="3.02" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="586.09" cy="232.04" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="4.03" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="255.79" r="3.83" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="156.84" r="6.14" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="556.92" cy="182.53" r="8.16" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="7.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="422.64" cy="552.77" r="5.94" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r=".7" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r=".5" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="8.06" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="1.71" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="2.42" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r=".6" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="5.74" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="9.57" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="9.97" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="579.96" cy="583.21" r="9.57" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="3.93" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="460.09" cy="571.94" r="9.87" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="1.91" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="429.32" cy="438.86" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="579.96" cy="583.21" r="9.97" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="6.44" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="7.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="579.96" cy="583.21" r="8.46" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="3.93" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="636.84" cy="479.86" r="8.96" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="9.97" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="630.62" cy="534.19" r="1.81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="1.41" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="636.84" cy="479.86" r="8.56" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="477.36" cy="394.63" r="4.73" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="569.19" cy="524.02" r="2.42" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="6.44" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="392.47" cy="119.19" r="7.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="4.43" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="9.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="1.51" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="392.47" cy="119.19" r="8.36" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="253.27" cy="156.84" r="2.22" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="284.53" cy="168.54" r="2.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="3.52" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="4.13" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="392.47" cy="119.19" r="2.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="348.35" cy="207.8" r="7.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="2.62" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="480.54" cy="152.84" r="7.85" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="392.47" cy="119.19" r="1.51" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="470.65" cy="207.24" r="2.82" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="480.54" cy="152.84" r="2.32" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="531.67" cy="156.84" r="8.16" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="392.47" cy="119.19" r=".5" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="480.54" cy="152.84" r=".81" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="7.65" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="8.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="4.13" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="1.01" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="4.53" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="369.41" cy="616.11" r="6.44" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="9.26" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="9.57" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="8.66" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="371.6" cy="650.86" r="7.75" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="7.55" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="6.54" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="8.86" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="432.99" cy="648.64" r="1.11" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="399.76" cy="595.51" r="3.12" fill=""none"" stroke=""none"" stroke-miterlimit="10"/><circle cx="385.72" cy="630.46" r="9.67" fill=""none"" stroke=""none"" stroke-miterlimit="10"/></g><g><g><circle cx="392.36" cy="20.47" r="20.47" fill=""none""/><rect x="386.23" y="33.45" width="12.25" height="49.04" fill=""none""/></g><g><circle cx="206.98" cy="69.31" r="20.47" fill=""none""/><rect x="219.61" y="77.26" width="12.25" height="49.04" transform="translate(-20.65 126.5) rotate(-30)" fill=""none""/></g><g><circle cx="70.86" cy="204.3" r="20.47" fill=""none""/><rect x="97.21" y="198.52" width="12.25" height="49.04" transform="translate(-141.5 201.01) rotate(-60)" fill=""none""/></g><g><circle cx="20.47" cy="389.26" r="20.47" fill=""none""/><rect x="51.84" y="364.74" width="12.25" height="49.04" transform="translate(-331.29 447.23) rotate(-90)" fill=""none""/></g><g><circle cx="69.31" cy="574.64" r="20.47" fill=""none""/><rect x="95.66" y="531.37" width="12.25" height="49.04" transform="translate(-328.73 921.98) rotate(-120)" fill=""none""/></g><g><circle cx="204.3" cy="710.76" r="20.47" fill=""none""/><rect x="216.92" y="653.76" width="12.25" height="49.04" transform="translate(77.07 1377.22) rotate(-150)" fill=""none""/></g><g><circle cx="389.26" cy="761.15" r="20.47" fill=""none""/><rect x="383.13" y="699.13" width="12.25" height="49.04" transform="translate(778.52 1447.3) rotate(-180)" fill=""none""/></g><g><circle cx="574.64" cy="712.31" r="20.47" fill=""none""/><rect x="549.76" y="655.31" width="12.25" height="49.04" transform="translate(1377.22 990.64) rotate(150)" fill=""none""/></g><g><circle cx="710.76" cy="577.32" r="20.47" fill=""none""/><rect x="672.16" y="534.05" width="12.25" height="49.04" transform="translate(1501.16 250.45) rotate(120)" fill=""none""/></g><g><circle cx="761.15" cy="392.36" r="20.47" fill=""none""/><rect x="717.53" y="367.84" width="12.25" height="49.04" transform="translate(1116.01 -331.29) rotate(90)" fill=""none""/></g><g><circle cx="712.31" cy="206.98" r="20.47" fill=""none""/><rect x="673.71" y="201.21" width="12.25" height="49.04" transform="translate(535.41 -475.89) rotate(60)" fill=""none""/></g><g><circle cx="577.32" cy="70.86" r="20.47" fill=""none""/><rect x="552.45" y="78.81" width="12.25" height="49.04" transform="translate(126.5 -265.44) rotate(30)" fill=""none""/></g></g></g></svg>
                        </div>
                        """
                    ),
                    class_ = "home-text"
                ),
                ui.div(
                    ui.HTML(
                            """
                            <h1 class="more-title">
                            <span class='span-more'>I have...</span>
                            </h1>
                            <div class="controlcontainer">
                                    <button class="button leftarrow">
                                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 640"><path d="M201.4 297.4C188.9 309.9 188.9 330.2 201.4 342.7L361.4 502.7C373.9 515.2 394.2 515.2 406.7 502.7C419.2 490.2 419.2 469.9 406.7 457.4L269.3 320L406.6 182.6C419.1 170.1 419.1 149.8 406.6 137.3C394.1 124.8 373.8 124.8 361.3 137.3L201.3 297.3z"></path></svg>
                                    </button>
                                    <p class="carousellabel" id="labeltext">Articles</p>
                                    <button class="button rightarrow">
                                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 640"><path d="M439.1 297.4C451.6 309.9 451.6 330.2 439.1 342.7L279.1 502.7C266.6 515.2 246.3 515.2 233.8 502.7C221.3 490.2 221.3 469.9 233.8 457.4L371.2 320L233.9 182.6C221.4 170.1 221.4 149.8 233.9 137.3C246.4 124.8 266.7 124.8 279.2 137.3L439.2 297.3z"></path></svg>
                                    </button>
                                    <script src="controls.js"></script>
                                </div>
                            <div class="carouselcontainer">
                                <input type="radio" name="position" checked />
                                <input type="radio" name="position" />
                                <input type="radio" name="position" />
                                <input type="radio" name="position" />
                                <div class="carousel">
                                <div class="item" page="arolit">
                                    <div class="header">
                                        <div class="img-box">
                                        <svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" class="arolitico" viewBox="0 0 640 640">                    <path d="M480 576L192 576C139 576 96 533 96 480L96 160C96 107 139 64 192 64L496 64C522.5 64 544 85.5 544 112L544 400C544 420.9 530.6 438.7 512 445.3L512 512C529.7 512 544 526.3 544 544C544 561.7 529.7 576 512 576L480 576zM192 448C174.3 448 160 462.3 160 480C160 497.7 174.3 512 192 512L448 512L448 448L192 448zM224 216C224 229.3 234.7 240 248 240L424 240C437.3 240 448 229.3 448 216C448 202.7 437.3 192 424 192L248 192C234.7 192 224 202.7 224 216zM248 288C234.7 288 224 298.7 224 312C224 325.3 234.7 336 248 336L424 336C437.3 336 448 325.3 448 312C448 298.7 437.3 288 424 288L248 288z"></path></svg>
                                        </div>
                                        <span class="title">AROLit</span>
                                    </div>
                                    <div class="content">
                                        <p>Query PubMed for articles of interest and then scrape potential primers from them.</p>
                                        <a class="action-button" id="arolitcard" href="#">Learn more</a>
                                    </div>
                                    </div>
                                    <div class="item" page="isop">
                                        <div class="header">
                                            <div class="img-box">
                                            <svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" class="isopico" viewBox="0 0 640 640">                    <path d="M295.4 37L310.2 73.8L347 88.6C350 89.8 352 92.8 352 96C352 99.2 350 102.2 347 103.4L310.2 118.2L295.4 155C294.2 158 291.2 160 288 160C284.8 160 281.8 158 280.6 155L265.8 118.2L229 103.4C226 102.2 224 99.2 224 96C224 92.8 226 89.8 229 88.6L265.8 73.8L280.6 37C281.8 34 284.8 32 288 32C291.2 32 294.2 34 295.4 37zM142.7 105.7L164.2 155.8L214.3 177.3C220.2 179.8 224 185.6 224 192C224 198.4 220.2 204.2 214.3 206.7L164.2 228.2L142.7 278.3C140.2 284.2 134.4 288 128 288C121.6 288 115.8 284.2 113.3 278.3L91.8 228.2L41.7 206.7C35.8 204.2 32 198.4 32 192C32 185.6 35.8 179.8 41.7 177.3L91.8 155.8L113.3 105.7C115.8 99.8 121.6 96 128 96C134.4 96 140.2 99.8 142.7 105.7zM496 368C502.4 368 508.2 371.8 510.7 377.7L532.2 427.8L582.3 449.3C588.2 451.8 592 457.6 592 464C592 470.4 588.2 476.2 582.3 478.7L532.2 500.2L510.7 550.3C508.2 556.2 502.4 560 496 560C489.6 560 483.8 556.2 481.3 550.3L459.8 500.2L409.7 478.7C403.8 476.2 400 470.4 400 464C400 457.6 403.8 451.8 409.7 449.3L459.8 427.8L481.3 377.7C483.8 371.8 489.6 368 496 368zM492 64C503 64 513.6 68.4 521.5 76.2L563.8 118.5C571.6 126.4 576 137 576 148C576 159 571.6 169.6 563.8 177.5L475.6 265.7L374.3 164.4L462.5 76.2C470.4 68.4 481 64 492 64zM76.2 462.5L340.4 198.3L441.7 299.6L177.5 563.8C169.6 571.6 159 576 148 576C137 576 126.4 571.6 118.5 563.8L76.2 521.5C68.4 513.6 64 503 64 492C64 481 68.4 470.4 76.2 462.5z"></path></svg><span class="title">iSOP</span>
                                            </div>
                                        </div>
                                        <div class="content">
                                            <p>
                                            Generate all <i>in silico primers</i> for a genome and score them based on a combination of important primer design parameters.
                                            </p>
                                            <a class="action-button" id="isopcard" href="#">Learn more</a>
                                        </div>
                                        </div>
                                    <div class="item" page="combos">
                                        <div class="header">
                                            <div class="img-box">
                                            <svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" class="combosico" viewBox="0 0 640 640">                    <path d="M96 160C96 124.7 124.7 96 160 96L480 96C515.3 96 544 124.7 544 160L544 480C544 515.3 515.3 544 480 544L160 544C124.7 544 96 515.3 96 480L96 160zM448 416C448 398.3 433.7 384 416 384C398.3 384 384 398.3 384 416C384 433.7 398.3 448 416 448C433.7 448 448 433.7 448 416zM224 256C241.7 256 256 241.7 256 224C256 206.3 241.7 192 224 192C206.3 192 192 206.3 192 224C192 241.7 206.3 256 224 256z"></path></svg><span page="combos" class="title">Primer Combinations</span>
                                            </div>
                                        </div>
                                        <div class="content">
                                            <p>Extract primer gene locations and generate the best primer combinations based on their scores.</p>
                                            <a class="action-button" id="comboscard" href="#">Learn more</a>
                                        </div>
                                        </div>
                                    <div class="item" page="utilities">
                                        <div class="header">
                                            <div class="img-box">
                                            <svg xmlns="http://www.w3.org/2000/svg" width="75" height="75" class="utilitiesico" viewBox="0 0 640 640">                    <path d="M240 120L240 160L400 160L400 120C400 115.6 396.4 112 392 112L248 112C243.6 112 240 115.6 240 120zM192 160L192 120C192 89.1 217.1 64 248 64L392 64C422.9 64 448 89.1 448 120L448 160L476.1 160C488.8 160 501 165.1 510 174.1L561.9 226C570.9 235 576 247.2 576 259.9L576 336L440 336L440 320C440 306.7 429.3 296 416 296C402.7 296 392 306.7 392 320L392 336L248 336L248 320C248 306.7 237.3 296 224 296C210.7 296 200 306.7 200 320L200 336L64 336L64 259.9C64 247.2 69.1 235 78.1 226L130 174.1C139 165.1 151.2 160 163.9 160L192 160zM64 480L64 384L200 384L200 400C200 413.3 210.7 424 224 424C237.3 424 248 413.3 248 400L248 384L392 384L392 400C392 413.3 402.7 424 416 424C429.3 424 440 413.3 440 400L440 384L576 384L576 480C576 515.3 547.3 544 512 544L128 544C92.7 544 64 515.3 64 480z"></path></svg>
                                            </div>
                                            <span class="title">Utilities</span>
                                        </div>
                                        <div class="content">
                                            <p>Various tools for data preprocessing and visualization.</p>
                                            <a class="action-button" id="utilitiescard" href="#">Learn more</a>
                                        </div>
                                        </div>
                                </div>
                        </div>"""
                ),
                class_ = "home-more"
            ),
            ui.div(

                class_ = "home-more"
            ),
        id = "home-page"
        ),
    ),

    # AROLit
    ui.nav_panel(
        "AROLit", 
        ui.div(
            ui.HTML("""<p><b><span class=span-arolit arolit-title>AROLit</span></b> - or <b>A</b>utomatic <b>R</b>etrieval of <b>O</b>ligonucleotides from <b>Lit</b>erature - is a tool for scraping PDFs for primers used in experiments, namely PCR. If the user has a set of search terms, AROLit can also query PubMed for articles matching those terms and download the associated metadata. Posteriorly, the retrieval of the full text PDFs can be done by uploading the metadata file to a reference manager such as Zotero.</p>"""),
            style="margin-bottom: 20px; text-align: justify;"
            ),
        ui.navset_bar(
            ui.nav_panel(
                     'Fetch articles',
                     ui.card(
                         ui.card_header("Fetch articles from PubMed based on search terms"),
                         ui.layout_columns(
                             ui.panel_well(
                                 ui.input_text("search_terms", "Enter search terms:", placeholder="e.g. 'Ebola' AND 'PCR'"),
                                 ui.input_action_button("fetch_articles", "Fetch Articles"),
                                 ui.output_text("nbib_status", inline=False, container=ui.tags.pre),
                             ),
                             ui.div(
                                 ui.h4("Articles Metadata", style="margin-bottom: 20px; font-weight: bold;"),
                                 ui.output_ui('metadata')
                             )
                         )
                     )

                 ),
            ui.nav_panel(
                'Primer scraping',
                ui.card(
                    ui.card_header('Scrape primers from literature'),
                    ui.layout_columns(
                        ui.panel_well(
                            ui.input_file('for_scraping','Upload PDF files for primer retrieval', multiple=True, accept=[".pdf"]),
                            ui.input_checkbox('sql_db', 'Save results to SQL database?'),
                            ui.input_action_button("scrape_primers", "Scrape Primers"),
                            ui.output_text("scrape_status", inline=False, container=ui.tags.pre),
                        ),
                        ui.div(
                            ui.h4("Scraped Primers", style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui("scraped_primers")
                        )
                    )
                )
                ),
            ui.nav_panel(
                'Primer validation',
                ui.card(
                    ui.card_header('Validate primers against the genome of interest'),
                    ui.layout_columns(
                        ui.panel_well(
                            ui.input_file('for_validation','Upload database file for primer validation',  accept=[".db"]),
                            ui.input_select(
                                "blast_method",
                                "Run BLAST against:",
                                choices={
                                    "genome_ref": "Reference genome",
                                    "genome_aligned": "Reference genome in an alignment"
                                },
                                selected="genome_ref"
                            ),
                            ui.panel_conditional(
                                "input.blast_method == 'genome_ref'",
                                ui.input_file("validate_ref", "Upload reference genome file:", accept=[".fasta", ".fa"])
                            ),
                            ui.panel_conditional(
                                "input.blast_method == 'genome_aligned'",
                                ui.input_file("validate_align", "Upload aligned reference genome file:", accept=[".fasta", ".fa"])
                            ),
                            ui.input_action_button("validate_primers", "Validate Primers"),                      
                            ui.output_text("validate_status", inline=False, container=ui.tags.pre),
                        ),
                        ui.div(
                            ui.h4("Valid Primers", style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui("valid_primers")
                        )
                    )
                )
                ),                
                 title=None,
        ),
        value='arolit'
    ),

    # iSOP
    ui.nav_panel(
        "iSOP",
        ui.div(
            ui.HTML("""<p><b><span class=span-isop isop-title>iSOP</span></b> - or <i><b>I</b>n <b>S</b>ilico</i> <b>O</b>ligonucleotides designed in <b>P</b>ython - allows users to generate every single possible primer sequence of a predetermined size for a given genome, even if the genome has been extracted from an alignment.</p>"""),
            style="margin-bottom: 20px; text-align: justify;"
        ),
        ui.navset_bar(
            ui.nav_panel(
                'Generate primers',
                ui.card(
                    ui.card_header(
                        ui.div(
                            'Generate ',
                            ui.tags.i('in silico'),
                            ' primers from a viral genome.'
                        )
                    ),
                    ui.layout_columns(
                        # Input column
                        ui.panel_well(
                            ui.input_file("genome_file", "Upload genome (FASTA)", accept=[".fasta", ".fa"]),
                            ui.input_text("primer_name", "Primer base name:", placeholder="e.g. Ebola"),
                            ui.input_numeric("min_length", "Minimum primer length:", 17, min=0, max=None),
                            ui.input_numeric("max_length", "Maximum primer length:", 18, min=0, max=None),
                            ui.input_checkbox("blast_filter", "Run BLAST for reference in an alignment?", value=False),
                            ui.panel_conditional(
                                "input.blast_filter == true",
                                ui.card(
                                    ui.input_file("align_ref", "Upload reference genome in the alignment", accept=[".fasta", ".fa"]),
                                )
                            ),
                            ui.input_action_button("generate_primers", "Generate Primers"),
                            ui.output_text("blast_status", inline=False, container=ui.tags.pre),
                        ),
                        # Output column
                        ui.div(
                            ui.h4(ui.tags.i("In silico"), " Primers", style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui("primer_outputs"),
                            class_='table-container'
                        ),
                    ),
                ),
            ),
            title=None,
            id='isop_tabs',
            navbar_options=ui.navbar_options(underline=True)
        ),
        value='isop'
    ),

    # Primer Combinations
    ui.nav_panel(
        "Primer Combinations",
        ui.div(
            ui.HTML("""<p><b><span class=span-combos combos-title>Primer Combinations</span></b> - a tool for scoring primers based on important parameters, such as GC content%, melting temperature, and more. In this tab, it's also possible to extract gene information from a sequence and map it to both a genome in an alignment and the primers themselves. Finally, the user can generate all possible valid primer combinations and evaluate their characteristics (Conservation scores & probability of <b>not</b> forming unwanted structures, such as hairpins, heterodimers or homodimers).</p>"""),
            style="margin-bottom: 20px; text-align: justify;"
        ),
        ui.navset_bar(
            ui.nav_panel(
                'Conservation Scores',
                ui.card(
                    ui.card_header("Calculate conservation scores for primers."),
                    ui.layout_columns(
                        ui.panel_well(
                            ui.input_file("align_file", "Upload alignment file of genomes (FASTA)", accept=[".fasta", ".fa"]),
                            ui.input_select(
                                "upload_method",
                                "Select primer dictionary origin:",
                                choices={
                                    "current": "Use the current generated primers",
                                    "upload": "Upload CSV file with primers"
                                },
                                selected="current"
                            ),
                            ui.panel_conditional(
                                "input.upload_method == 'upload'",
                                ui.input_file("primer_csv", "Upload primer CSV file:", accept=[".csv"])
                            ),
                            ui.input_action_button("calculate_CS", "Calculate Conservation Scores"),
                            ui.output_text("CS_status", inline=False, container=ui.tags.pre),
                        ),
                        # Output column
                        ui.div(
                            ui.h4("Conservation Scores",
                                  style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui("CS_output"),
                        )
                    )
                ),
                value="conservation_scores"
            ),            
            ui.nav_panel(
                'Gene extraction',
                ui.card(
                    ui.card_header(
                        ui.div(
                            'Extract genes from a reference genome and map the loci to a CSV containing primers',
                        )
                    ),
                    ui.layout_columns(
                        ui.panel_well(
                            ui.input_select(
                                "retrieve_method",
                                "Select your gene data retrieval method:",
                                choices={
                                    "auto": "Automatically retrieve from NCBI",
                                    "upload": "Upload GenBank file"
                                },
                                selected="auto"
                            ),
                            ui.panel_conditional(
                                "input.retrieve_method == 'upload'",
                                ui.input_file("genbank_file", "Upload GenBank file:", accept=[".gb", ".gbk"])
                            ),
                            ui.panel_conditional(
                                "input.retrieve_method == 'auto'",
                                ui.input_text("ncbi_id", "Enter NCBI ID:"),
                                ui.input_text('email', 'Enter your email address:')
                            ),
                            ui.input_checkbox("map_to", "Map the loci to the reference in an alignment?", value=False),
                            ui.panel_conditional(
                                "input.map_to == true",
                                ui.card(
                                    ui.input_file("reference_file", "Upload reference sequence:", accept=[".fasta", ".fa"]),
                                    ui.input_file("alignment_file", "Upload alignment file:", accept=[".fasta", ".fa"]),
                                )
                            ),
                            ui.input_action_button("retrieve_loci", "Retrieve gene loci"),
                            ui.output_text("loci_status", inline=False, container=ui.tags.pre),
                        ),
                        # Output column
                        ui.div(
                            ui.h4("Loci Information", style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui("lociref_output"),
                            ui.output_ui("locialign_output")
                        )
                    ),
                )
            ),
            ui.nav_panel(
                'Annotate primers',
                ui.card(
                    ui.card_header('Annotate primers with loci information'),
                    ui.layout_columns(
                        # Input column
                        ui.panel_well(
                            ui.input_file('loci_file', 'Upload file with loci information:', accept=[".csv"]),
                            ui.input_file('primer_file', 'Upload file with primers to annotate', accept=[".csv"]),
                            ui.help_text(
                                ui.tags.i(class_="fa-solid fa-triangle-exclamation"),
                                ' Please make sure the primers file has the same structure as the CSV obtained in the ',
                                ui.input_action_link('cs_switch', 'Conservation Scores'),         
                                ' section.',
                                ui.tags.br(),
                            ),
                            ui.input_action_button('annotate_primers', 'Annotate Primers'),
                            ui.output_text("annotated_status", inline=False, container=ui.tags.pre),
                        ),
                        # Output column
                        ui.div(
                            ui.h4('Annotated Primers', style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui('annotated_primers')
                        )
                    )
                ),
                value="annotate_primers"
            ),
            ui.nav_panel(
                'Primer combinations',
                ui.card(
                    ui.card_header('Explore primer combinations'),
                    ui.layout_columns(
                        # Input column
                        ui.panel_well(
                            ui.input_file('primers_ready', 'Upload file with primers:', accept=[".csv"]),
                            ui.help_text(
                                ui.tags.i(class_="fa-solid fa-triangle-exclamation"),
                                ' Please make sure the primers file has the same structure as the output from the ',
                                ui.input_action_link('annotation_switch', 'Annotate primers'),         
                                ' section.',
                                ui.tags.br(),
                            ),                            
                            ui.input_switch('custom_parameters', 'Customize combination parameters'),
                            ui.panel_conditional(
                                "input.custom_parameters == true",
                                ui.card(
                                    ui.card_header('Custom Parameters'),
                                    ui.layout_columns(
                                    ui.input_numeric('amp_min', 'Minimum Amplicon Length', value=100, min=1),
                                    ui.input_numeric('amp_max', 'Maximum Amplicon Length', value=1000, min=1),
                                    ),
                                    ui.layout_columns(
                                    ui.input_numeric('gc_min', 'Minimum GC Content (%)', value=40, min=0, max=100),
                                    ui.input_numeric('gc_max', 'Maximum GC Content (%)', value=60, min=0, max=100),
                                    ),
                                    ui.layout_columns(
                                    ui.input_numeric('mt_min', 'Minimum Melting Temperature (°C)', value=55, min=0),
                                    ui.input_numeric('mt_max', 'Maximum Melting Temperature (°C)', value=65, min=0),
                                    ),
                                    ui.input_numeric('min_cs', 'Conservation Score threshold (%)', value=99, min=0, max=100),
                                )
                            ),
                            ui.input_checkbox('show_histograms', 'Show histograms for GC content and Melting Temperature distribution'),
                            ui.input_action_button('explore_combinations', 'Calculate primer combinations'),
                            ui.output_text("combinations_status", inline=False, container=ui.tags.pre),
                        ),
                        # Output column
                        ui.div(
                            ui.h4('Primer Combinations', style="margin-bottom: 20px; font-weight: bold;"),
                            ui.output_ui('primer_combinations')
                        )
                    )
                ),
                ui.panel_conditional(
                    "input.show_histograms == true",
                    ui.card(
                        ui.card_header('Histograms'),
                        ui.layout_columns(
                            ui.output_plot('gc_histogram'),
                            ui.output_plot('tm_histogram')
                        )
                    )
                )
            ),
            title=None,
            id='combinations_tabs'
        ),
        value='primer_combinations'
    ),

    # Utilities
    ui.nav_panel(
        "Utilities", 
        ui.div(
            ui.HTML("""<p><b><span class=span-utilities utilities-title>Utilities</span></b> - assortment of tools to help pre-process raw data or prepare inputs for the ViruScope tools.</p>"""),
            style="margin-bottom: 20px; text-align: justify;"
        ),
        ui.navset_pill_list(
            ui.nav_control(
                ui.tags.div(
                ui.tags.strong("Tools"),
                class_="nav-label"
                )
            ),
            ui.nav_panel("Fetch DOIs and PMIDs",
                         ui.card(
                           ui.card_header("Compile a list of all DOIs/PMIDs from a .nbib file"),
                           ui.card_body(
                               ui.layout_columns(
                                   ui.panel_well(
                                       ui.input_file("nbib_file", "Upload .nbib file:", accept=[".nbib"]),
                                       ui.input_action_button("fetch_dois", "Fetch DOIs/PMIDs", width='fit-content'),
                                       ui.output_text("fetch_status", inline=False, container=ui.tags.pre)
                                    ),
                                      ui.div(
                                        ui.h4("List of DOIs",style="margin-bottom: 20px; font-weight: bold;"),
                                        ui.output_ui("doi_output"),
                                        style='border-right: 1px solid #ccc;'
                                      ),
                                      ui.div(
                                        ui.h4("List of PMIDs",style="margin-bottom: 20px; font-weight: bold;"),
                                        ui.output_ui("pmid_output"),
                                      )
                               )
                           )
                         ),
                         ),
            ui.nav_panel("Count sequences",
                         ui.card(
                           ui.card_header("Count sequences in a FASTA file"),
                           ui.card_body(
                             ui.input_file("fasta_file", "Upload FASTA file", accept=[".fasta", ".fa"]),
                             ui.input_action_button("count_sequences", "Count Sequences", width='fit-content'),
                             ui.tags.div(
                               ui.tags.b("Total sequences in the file: "),
                               ui.output_text("sequence_count", inline=True),
                               style="display: inline-block;"
                             )
                           )
                         ),
                         ),
            ui.nav_panel("Reorder sequences", 
                        ui.card(
                           ui.card_header("Reorder sequences in a FASTA file"),
                           ui.card_body(
                             ui.input_file("fasta_file_reorder", "Upload FASTA file:", accept=[".fasta", ".fa"]),
                             ui.input_text("reference_id", "Reference ID:", placeholder="e.g. NC_002549.1"),
                             ui.output_ui("download_reordered_ui"),
                           )
                         ),
                         ),
            ui.nav_panel("Remove duplicates", 
                        ui.card(
                           ui.card_header("Remove duplicate sequences from a FASTA file"),
                           ui.card_body(
                             ui.input_file("fasta_file_dupes", "Upload FASTA file:", accept=[".fasta", ".fa"]),
                            ui.help_text(
                                ui.tags.i(class_="fa-solid fa-triangle-exclamation"),
                                ' Reference sequences should be the first entry in the FASTA file to avoid getting deleted.',
                                ui.tags.br(),
                            ),
                             ui.output_ui("download_nodupes_ui"),
                           )
                         ),
                         ),
            ui.nav_panel("Extract sequence",
                        ui.card(
                        ui.card_header("Extract a specific sequence from a FASTA file"),
                        ui.card_body(
                             ui.input_file("fasta_file_extract", "Upload FASTA file:", accept=[".fasta", ".fa"]),
                             ui.input_text("sequence_id", "Sequence ID:", placeholder="e.g. NC_002549.1"),
                             ui.output_ui("download_extract_ui")
                           )
                         ),
                         ),
            ui.nav_panel("Venn diagram",
                        ui.card(
                        ui.card_header("Generate a Venn diagram of primer sequences"),
                        ui.layout_columns(
                            ui.div(
                                ui.input_file("fasta_file_venn1", "Upload FASTA file 1:", accept=[".fasta", ".fa"]),
                                ui.input_file("fasta_file_venn2", "Upload FASTA file 2:", accept=[".fasta", ".fa"]),
                                ui.input_action_button("generate_venn", "Generate Venn Diagram", width='fit-content'),
                                style='border-right: 1px solid #ccc;',
                            ),
                            ui.output_plot("venn_diagram"),
                            col_widths=(4,8),
                           ),
                         ),
                         ),
            widths=(2,10)
        ),
        value='utilities'

    ),

    # Dark mode
    ui.nav_control(ui.input_dark_mode(id="mode", mode="light")),

    # Logo
    title=ui.div(
        ui.tags.a(
            ui.output_ui("toggle_logo")
        )
    ),
    footer=ui.tags.footer(
        ui.tags.ul(
            {"class": "socials"},
            ui.tags.li(
                {"class": "icon-content"},
                ui.tags.a(
                    {"href": "https://github.com/anasfplima", 'aria-label':'Github', "data-social": "github"},
                    ui.tags.div(class_ = "filled"),
                    ui.tags.svg(
                        {'xmlns':"http://www.w3.org/2000/svg",
                        'width':"16",
                        'height':"16",
                        'fill':"currentColor",
                        'viewBox':"0 0 640 640"},
                        ui.tags.Tag('path', d="M237.9 461.4C237.9 463.4 235.6 465 232.7 465C229.4 465.3 227.1 463.7 227.1 461.4C227.1 459.4 229.4 457.8 232.3 457.8C235.3 457.5 237.9 459.1 237.9 461.4zM206.8 456.9C206.1 458.9 208.1 461.2 211.1 461.8C213.7 462.8 216.7 461.8 217.3 459.8C217.9 457.8 216 455.5 213 454.6C210.4 453.9 207.5 454.9 206.8 456.9zM251 455.2C248.1 455.9 246.1 457.8 246.4 460.1C246.7 462.1 249.3 463.4 252.3 462.7C255.2 462 257.2 460.1 256.9 458.1C256.6 456.2 253.9 454.9 251 455.2zM316.8 72C178.1 72 72 177.3 72 316C72 426.9 141.8 521.8 241.5 555.2C254.3 557.5 258.8 549.6 258.8 543.1C258.8 536.9 258.5 502.7 258.5 481.7C258.5 481.7 188.5 496.7 173.8 451.9C173.8 451.9 162.4 422.8 146 415.3C146 415.3 123.1 399.6 147.6 399.9C147.6 399.9 172.5 401.9 186.2 425.7C208.1 464.3 244.8 453.2 259.1 446.6C261.4 430.6 267.9 419.5 275.1 412.9C219.2 406.7 162.8 398.6 162.8 302.4C162.8 274.9 170.4 261.1 186.4 243.5C183.8 237 175.3 210.2 189 175.6C209.9 169.1 258 202.6 258 202.6C278 197 299.5 194.1 320.8 194.1C342.1 194.1 363.6 197 383.6 202.6C383.6 202.6 431.7 169 452.6 175.6C466.3 210.3 457.8 237 455.2 243.5C471.2 261.2 481 275 481 302.4C481 398.9 422.1 406.6 366.2 412.9C375.4 420.8 383.2 435.8 383.2 459.3C383.2 493 382.9 534.7 382.9 542.9C382.9 549.4 387.5 557.3 400.2 555C500.2 521.8 568 426.9 568 316C568 177.3 455.5 72 316.8 72zM169.2 416.9C167.9 417.9 168.2 420.2 169.9 422.1C171.5 423.7 173.8 424.4 175.1 423.1C176.4 422.1 176.1 419.8 174.4 417.9C172.8 416.3 170.5 415.6 169.2 416.9zM158.4 408.8C157.7 410.1 158.7 411.7 160.7 412.7C162.3 413.7 164.3 413.4 165 412C165.7 410.7 164.7 409.1 162.7 408.1C160.7 407.5 159.1 407.8 158.4 408.8zM190.8 444.4C189.2 445.7 189.8 448.7 192.1 450.6C194.4 452.9 197.3 453.2 198.6 451.6C199.9 450.3 199.3 447.3 197.3 445.4C195.1 443.1 192.1 442.8 190.8 444.4zM179.4 429.7C177.8 430.7 177.8 433.3 179.4 435.6C181 437.9 183.7 438.9 185 437.9C186.6 436.6 186.6 434 185 431.7C183.6 429.4 181 428.4 179.4 429.7z"),
                    ),
                    ui.tags.span("GitHub", class_="tooltip")
                 ),               
            ),
            ui.tags.li(
                {"class": "icon-content"},
                ui.tags.a(
                    {"href": "https://viruscope.jc-biotechaiteam.com/wp/", 'aria-label':'Github', "data-social": "db"},
                    ui.tags.div(class_ = "filled"),
                    ui.tags.svg(
                        {'xmlns':"http://www.w3.org/2000/svg",
                        'width':"16",
                        'height':"16",
                        'fill':"currentColor",
                        'viewBox':"80 25 366 366"},
                        ui.tags.Tag('path', d="M259.06 36.96h5.74v23h-5.74zm11.182-10.888a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm72.893 32.254 4.971 2.87-11.5 19.919-4.97-2.87zm13.841-8.625a9.6 9.6 0 0 1-3.514 13.113 9.6 9.6 0 0 1-13.114-3.513 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm48.283 69.188 2.87 4.97-19.919 11.5-2.87-4.97zm15.014-5.355a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm8.51 88.846v5.74h-23v-5.74zm14.39-1.92a9.6 9.6 0 0 1-3.513 13.113 9.6 9.6 0 0 1-13.114-3.513 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm-35.765 86.004-2.87 4.97-19.919-11.5 2.87-4.97zm12.137.73a9.6 9.6 0 0 1-3.513 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.113 3.514zm-72.693 61.392-4.97 2.87-11.5-19.919 4.97-2.87zm8.86 1.905a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm-92.361 21.62h-5.74v-23h5.74zm5.435 1.28a9.6 9.6 0 0 1-3.513 13.115 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.113 3.514zM179.27 350.73l-4.971-2.87 11.5-19.919 4.97 2.87zm2.782-.966a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.113 9.6 9.6 0 0 1 13.114 3.513zm-64.897-59.584-2.87-4.971 19.919-11.5 2.87 4.971zm1.6-4.25a9.6 9.6 0 0 1-3.514 13.115 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm-25.133-79.248v-5.74h23v5.74zm2.232-7.677a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.113-3.514 9.6 9.6 0 0 1 3.513-13.113 9.6 9.6 0 0 1 13.114 3.513zm19.151-76.405 2.87-4.971 19.919 11.5-2.87 4.97zm4.478-10.329a9.6 9.6 0 0 1-3.514 13.114 9.6 9.6 0 0 1-13.114-3.514 9.6 9.6 0 0 1 3.514-13.114 9.6 9.6 0 0 1 13.114 3.514zm56.074-51.797 4.97-2.87 11.5 19.919-4.97 2.87zm7.76-11.5a9.6 9.6 0 0 1-3.515 13.114 9.6 9.6 0 0 1-13.113-3.514 9.6 9.6 0 0 1 3.513-13.114 9.6 9.6 0 0 1 13.114 3.514zm135.935 251.808a4.53 4.53 0 0 1-1.658 6.188 4.53 4.53 0 0 1-6.188-1.658 4.53 4.53 0 0 1 1.658-6.188 4.53 4.53 0 0 1 6.188 1.658zm-5.155-15.95a1.46 1.46 0 0 1-.534 1.995 1.46 1.46 0 0 1-1.995-.534 1.46 1.46 0 0 1 .535-1.995 1.46 1.46 0 0 1 1.994.535zm25.134 14.253a.52.52 0 0 1-.19.71.52.52 0 0 1-.711-.19.52.52 0 0 1 .19-.71.52.52 0 0 1 .71.19zm-20.3 1.882a4.16 4.16 0 0 1-1.522 5.683 4.16 4.16 0 0 1-5.683-1.523 4.16 4.16 0 0 1 1.523-5.683 4.16 4.16 0 0 1 5.683 1.523zm-.943.545a3.07 3.07 45 0 1-1.124 4.194 3.07 3.07 45 0 1-4.194-1.124 3.07 3.07 45 0 1 1.124-4.194 3.07 3.07 45 0 1 4.194 1.124zm23.858-3.937a3.54 3.54 0 0 1-1.296 4.836 3.54 3.54 0 0 1-4.836-1.296 3.54 3.54 0 0 1 1.296-4.835 3.54 3.54 0 0 1 4.836 1.295zm-24.321 15.255a3.64 3.64 0 0 1-1.333 4.972 3.64 3.64 0 0 1-4.972-1.332 3.64 3.64 0 0 1 1.332-4.973 3.64 3.64 0 0 1 4.973 1.333zm1.32-11.813a4.06 4.06 0 0 1-1.486 5.546 4.06 4.06 0 0 1-5.546-1.486 4.06 4.06 0 0 1 1.486-5.546 4.06 4.06 0 0 1 5.546 1.486zm-.584 11.388a4.49 4.49 0 0 1-1.644 6.133 4.49 4.49 0 0 1-6.133-1.643 4.49 4.49 0 0 1 1.643-6.134 4.49 4.49 0 0 1 6.134 1.644zm.827-11.528a4.34 4.34 0 0 1-1.589 5.928 4.34 4.34 0 0 1-5.928-1.588 4.34 4.34 0 0 1 1.588-5.929 4.34 4.34 0 0 1 5.929 1.589zm-11.134-1.343a3.02 3.02 0 0 1-1.105 4.125 3.02 3.02 0 0 1-4.125-1.105 3.02 3.02 0 0 1 1.105-4.126 3.02 3.02 0 0 1 4.125 1.106zm6.714-15.031a2.12 2.12 45 0 1-.775 2.896 2.12 2.12 45 0 1-2.896-.776 2.12 2.12 45 0 1 .776-2.896 2.12 2.12 45 0 1 2.895.776zm1.068 18.309a.47.47 45 0 1-.172.642.47.47 45 0 1-.642-.172.47.47 45 0 1 .172-.642.47.47 45 0 1 .642.172zm.316 10.868a1.94 1.94 45 0 1-.71 2.65 1.94 1.94 45 0 1-2.65-.71 1.94 1.94 45 0 1 .71-2.65 1.94 1.94 45 0 1 2.65.71zm-7.647-14.406a3.54 3.54 0 0 1-1.296 4.835 3.54 3.54 0 0 1-4.836-1.295 3.54 3.54 0 0 1 1.296-4.836 3.54 3.54 0 0 1 4.836 1.296zm10.276 1.838a3.87 3.87 0 0 1-1.417 5.286 3.87 3.87 0 0 1-5.286-1.416 3.87 3.87 0 0 1 1.416-5.287 3.87 3.87 0 0 1 5.287 1.417zm-1.2 11.743a3.59 3.59 0 0 1-1.314 4.904 3.59 3.59 0 0 1-4.904-1.314 3.59 3.59 0 0 1 1.314-4.904 3.59 3.59 0 0 1 4.904 1.314zM242.177 86.662a.38.38 0 0 1-.14.52.38.38 0 0 1-.519-.14.38.38 0 0 1 .14-.519.38.38 0 0 1 .519.14zm-43.778 7.054a.24.24 0 0 1-.088.328.24.24 0 0 1-.328-.088.24.24 0 0 1 .088-.328.24.24 0 0 1 .328.088zm68.464-19.136a3.82 3.82 0 0 1-1.398 5.219 3.82 3.82 0 0 1-5.218-1.399 3.82 3.82 0 0 1 1.398-5.218 3.82 3.82 0 0 1 5.218 1.398zm-24.071 11.727a1.09 1.09 0 0 1-.4 1.49 1.09 1.09 0 0 1-1.488-.4 1.09 1.09 0 0 1 .399-1.489 1.09 1.09 0 0 1 1.489.4zm8.935 24.298a1.32 1.32 0 0 1-.483 1.803 1.32 1.32 0 0 1-1.803-.483 1.32 1.32 0 0 1 .483-1.804 1.32 1.32 0 0 1 1.803.484zM198.806 93.48a.71.71 0 0 1-.26.97.71.71 0 0 1-.97-.26.71.71 0 0 1 .26-.97.71.71 0 0 1 .97.26zm46.229-8.469a3.68 3.68 0 0 1-1.347 5.027 3.68 3.68 0 0 1-5.027-1.347 3.68 3.68 0 0 1 1.347-5.027 3.68 3.68 0 0 1 5.027 1.347zm6.614 25.638a1.23 1.23 0 0 1-.45 1.68 1.23 1.23 0 0 1-1.68-.45 1.23 1.23 0 0 1 .45-1.68 1.23 1.23 0 0 1 1.68.45zM204.2 138.345a3.64 3.64 0 0 1-1.332 4.972 3.64 3.64 0 0 1-4.973-1.332 3.64 3.64 0 0 1 1.333-4.973 3.64 3.64 0 0 1 4.972 1.333zm-4.944-45.124a1.23 1.23 0 0 1-.45 1.68 1.23 1.23 0 0 1-1.68-.45 1.23 1.23 0 0 1 .45-1.68 1.23 1.23 0 0 1 1.68.45zm53.008 17.074a1.94 1.94 45 0 1-.71 2.65 1.94 1.94 45 0 1-2.65-.71 1.94 1.94 45 0 1 .71-2.65 1.94 1.94 45 0 1 2.65.71zm-49.787 29.045a1.65 1.65 0 0 1-.604 2.254 1.65 1.65 0 0 1-2.254-.604 1.65 1.65 0 0 1 .604-2.254 1.65 1.65 0 0 1 2.254.604zm-35.732-.629a.94.94 45 0 1-.344 1.284.94.94 45 0 1-1.284-.344.94.94 45 0 1 .344-1.284.94.94 45 0 1 1.284.344zm-15.358 2.54a1.04 1.04 0 0 1-.38 1.42 1.04 1.04 0 0 1-1.421-.38 1.04 1.04 0 0 1 .38-1.421 1.04 1.04 0 0 1 1.421.38zm50.198-49.375a3.92 3.92 45 0 1-1.434 5.355 3.92 3.92 45 0 1-5.355-1.435 3.92 3.92 45 0 1 1.435-5.354 3.92 3.92 45 0 1 5.355 1.434zm-35.039 46.95a.71.71 0 0 1-.26.97.71.71 0 0 1-.97-.26.71.71 0 0 1 .26-.97.71.71 0 0 1 .97.26zm38.303-.856a4.39 4.39 0 0 1-1.606 5.997 4.39 4.39 0 0 1-5.997-1.607 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.996 1.607zm-37.117.171a2.08 2.08 0 0 1-.76 2.842 2.08 2.08 0 0 1-2.842-.762 2.08 2.08 0 0 1 .761-2.841 2.08 2.08 0 0 1 2.841.761zm33.611-46.125a3.64 3.64 0 0 1-1.332 4.973 3.64 3.64 0 0 1-4.973-1.333 3.64 3.64 0 0 1 1.333-4.972 3.64 3.64 0 0 1 4.972 1.332zm2.32 46.639a3.02 3.02 0 0 1-1.105 4.125 3.02 3.02 0 0 1-4.126-1.105 3.02 3.02 0 0 1 1.106-4.126 3.02 3.02 0 0 1 4.125 1.106zm162.194 77.588a1.13 1.13 0 0 1-.414 1.543 1.13 1.13 0 0 1-1.543-.413 1.13 1.13 0 0 1 .413-1.544 1.13 1.13 0 0 1 1.544.414zm-66.682-31.557a2.22 2.22 0 0 1-.813 3.032 2.22 2.22 0 0 1-3.032-.812 2.22 2.22 0 0 1 .812-3.033 2.22 2.22 0 0 1 3.033.813zm86.291-3.679a4.01 4.01 0 0 1-1.467 5.478 4.01 4.01 0 0 1-5.478-1.468 4.01 4.01 0 0 1 1.468-5.478 4.01 4.01 0 0 1 5.477 1.468zm-20.016 35.47a.66.66 0 0 1-.242.902.66.66 0 0 1-.901-.241.66.66 0 0 1 .241-.902.66.66 0 0 1 .902.242zm27.49-10.363a.85.85 0 0 1-.31 1.16.85.85 0 0 1-1.161-.31.85.85 0 0 1 .31-1.162.85.85 0 0 1 1.162.312zm-24.017 8.359a4.67 4.67 0 0 1-1.71 6.379 4.67 4.67 0 0 1-6.38-1.71 4.67 4.67 0 0 1 1.71-6.379 4.67 4.67 0 0 1 6.38 1.71zm16.708-33.56a4.2 4.2 0 0 1-1.537 5.737 4.2 4.2 0 0 1-5.738-1.538 4.2 4.2 0 0 1 1.538-5.737 4.2 4.2 0 0 1 5.737 1.537zm8.167 24.706a1.84 1.84 0 0 1-.673 2.513 1.84 1.84 0 0 1-2.514-.673 1.84 1.84 0 0 1 .674-2.514 1.84 1.84 0 0 1 2.513.674zm-7.233 30.71a3.97 3.97 0 0 1-1.453 5.424 3.97 3.97 0 0 1-5.423-1.453 3.97 3.97 0 0 1 1.453-5.423 3.97 3.97 0 0 1 5.423 1.453zm-18.578-21.316a3.59 3.59 0 0 1-1.314 4.904 3.59 3.59 0 0 1-4.904-1.314 3.59 3.59 0 0 1 1.314-4.904 3.59 3.59 0 0 1 4.904 1.314zm26.833-9.984a3.02 3.02 0 0 1-1.105 4.125 3.02 3.02 0 0 1-4.126-1.105 3.02 3.02 0 0 1 1.106-4.126 3.02 3.02 0 0 1 4.125 1.106zm-7.649 30.95a4.67 4.67 0 0 1-1.71 6.38 4.67 4.67 0 0 1-6.379-1.71 4.67 4.67 0 0 1 1.71-6.379 4.67 4.67 0 0 1 6.379 1.71zm-95.142-23.23a4.53 4.53 0 0 1-1.658 6.188 4.53 4.53 0 0 1-6.188-1.658 4.53 4.53 0 0 1 1.658-6.189 4.53 4.53 0 0 1 6.188 1.659zm8.289-28.723a3.54 3.54 0 0 1-1.296 4.836 3.54 3.54 0 0 1-4.836-1.296 3.54 3.54 0 0 1 1.296-4.836 3.54 3.54 0 0 1 4.836 1.296zm65.34 32.332a.9.9 0 0 1-.33 1.229.9.9 0 0 1-1.23-.33.9.9 0 0 1 .33-1.229.9.9 0 0 1 1.23.33zm-73.75-3.54a4.39 4.39 0 0 1-1.607 5.997 4.39 4.39 0 0 1-5.997-1.606 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.996 1.607zm43.91 46.714a4.63 4.63 0 0 1-1.696 6.325 4.63 4.63 0 0 1-6.324-1.695 4.63 4.63 0 0 1 1.694-6.324 4.63 4.63 0 0 1 6.325 1.694zm30.654-43.644a1.84 1.84 0 0 1-.674 2.513 1.84 1.84 0 0 1-2.513-.673 1.84 1.84 0 0 1 .673-2.514 1.84 1.84 0 0 1 2.514.674zm20.543 20.182a4.49 4.49 0 0 1-1.643 6.133 4.49 4.49 0 0 1-6.134-1.643 4.49 4.49 0 0 1 1.644-6.134 4.49 4.49 0 0 1 6.133 1.644zm-51.77 23.792a3.97 3.97 0 0 1-1.453 5.423 3.97 3.97 0 0 1-5.423-1.453 3.97 3.97 0 0 1 1.453-5.423 3.97 3.97 0 0 1 5.423 1.453zm.607-.35a4.67 4.67 0 0 1-1.71 6.38 4.67 4.67 0 0 1-6.379-1.71 4.67 4.67 0 0 1 1.71-6.38 4.67 4.67 0 0 1 6.379 1.71zm-43.858-46.743a4.49 4.49 0 0 1-1.643 6.133 4.49 4.49 0 0 1-6.134-1.643 4.49 4.49 0 0 1 1.644-6.134 4.49 4.49 0 0 1 6.133 1.644zm75.214 2.694a2.69 2.69 45 0 1-.985 3.674 2.69 2.69 45 0 1-3.674-.984 2.69 2.69 45 0 1 .984-3.675 2.69 2.69 45 0 1 3.675.985zm-35.158 46.244a.28.28 0 0 1-.103.383.28.28 0 0 1-.382-.103.28.28 0 0 1 .102-.382.28.28 0 0 1 .383.102zm7.709 37.072a1.13 1.13 0 0 1-.414 1.544 1.13 1.13 0 0 1-1.543-.414 1.13 1.13 0 0 1 .413-1.543 1.13 1.13 0 0 1 1.544.413zm-26.234-13.617a.8.8 0 0 1-.292 1.093.8.8 0 0 1-1.093-.293.8.8 0 0 1 .293-1.093.8.8 0 0 1 1.092.293zm21.556-25.205a3.78 3.78 0 0 1-1.384 5.164 3.78 3.78 0 0 1-5.163-1.384 3.78 3.78 0 0 1 1.383-5.163 3.78 3.78 0 0 1 5.164 1.383zm3.907 39.267a.24.24 0 0 1-.088.328.24.24 0 0 1-.328-.088.24.24 0 0 1 .088-.327.24.24 0 0 1 .328.087zm-25.87-13.827a.33.33 0 0 1-.12.45.33.33 0 0 1-.451-.12.33.33 0 0 1 .12-.451.33.33 0 0 1 .451.12zm1.403-23.95a2.79 2.79 0 0 1-1.021 3.81 2.79 2.79 0 0 1-3.811-1.02 2.79 2.79 0 0 1 1.02-3.812 2.79 2.79 0 0 1 3.812 1.022zm20.473-1.44a3.68 3.68 0 0 1-1.347 5.027 3.68 3.68 0 0 1-5.027-1.347 3.68 3.68 0 0 1 1.347-5.027 3.68 3.68 0 0 1 5.027 1.347zm-19.096 23.785a3.54 3.54 0 0 1-1.295 4.835 3.54 3.54 0 0 1-4.836-1.295 3.54 3.54 0 0 1 1.296-4.836 3.54 3.54 0 0 1 4.835 1.296zm-32.77-204.7a3.82 3.82 0 0 1-1.398 5.219 3.82 3.82 0 0 1-5.218-1.399 3.82 3.82 0 0 1 1.398-5.218 3.82 3.82 0 0 1 5.218 1.398zM266.05 75.05a2.88 2.88 45 0 1-1.054 3.935 2.88 2.88 45 0 1-3.934-1.055 2.88 2.88 45 0 1 1.054-3.934 2.88 2.88 45 0 1 3.934 1.054zm62.43 17.534a1.79 1.79 0 0 1-.654 2.445 1.79 1.79 0 0 1-2.445-.655 1.79 1.79 0 0 1 .655-2.445 1.79 1.79 0 0 1 2.445.655zm-47.022-12.527a1.89 1.89 0 0 1-.692 2.582 1.89 1.89 0 0 1-2.582-.692 1.89 1.89 0 0 1 .692-2.581 1.89 1.89 0 0 1 2.582.691zm22.634 13.744a.94.94 45 0 1-.344 1.284.94.94 45 0 1-1.284-.344.94.94 45 0 1 .344-1.284.94.94 45 0 1 1.284.344zM281.05 80.292a1.42 1.42 0 0 1-.52 1.94 1.42 1.42 0 0 1-1.94-.52 1.42 1.42 0 0 1 .52-1.94 1.42 1.42 0 0 1 1.94.52zm47.43 12.292a1.79 1.79 0 0 1-.655 2.445 1.79 1.79 0 0 1-2.445-.655 1.79 1.79 0 0 1 .655-2.445 1.79 1.79 0 0 1 2.445.655zm-22.137-.083a3.54 3.54 0 0 1-1.296 4.836 3.54 3.54 0 0 1-4.836-1.296 3.54 3.54 0 0 1 1.296-4.835 3.54 3.54 0 0 1 4.836 1.295zm-33.447-3.171a4.53 4.53 0 0 1-1.658 6.188 4.53 4.53 0 0 1-6.188-1.658 4.53 4.53 0 0 1 1.658-6.188 4.53 4.53 0 0 1 6.188 1.658zm7.66-8.753a.85.85 0 0 1-.311 1.161.85.85 0 0 1-1.161-.31.85.85 0 0 1 .31-1.162.85.85 0 0 1 1.162.311zm24.029 12.94a1.51 1.51 0 0 1-.553 2.062 1.51 1.51 0 0 1-2.063-.553 1.51 1.51 0 0 1 .553-2.062 1.51 1.51 0 0 1 2.063.552zm-33.897-2.912a1.98 1.98 0 0 1-.725 2.705 1.98 1.98 0 0 1-2.705-.725 1.98 1.98 0 0 1 .725-2.705 1.98 1.98 0 0 1 2.705.725zm-2.86-4.492a4.39 4.39 0 0 1-1.606 5.996 4.39 4.39 0 0 1-5.997-1.606 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.997 1.607zm-1.614-11.158a3.07 3.07 45 0 1-1.124 4.194 3.07 3.07 45 0 1-4.194-1.124 3.07 3.07 45 0 1 1.124-4.193 3.07 3.07 45 0 1 4.194 1.123zm16.914 4.137a3.82 3.82 0 0 1-1.398 5.219 3.82 3.82 0 0 1-5.218-1.399 3.82 3.82 0 0 1 1.398-5.218 3.82 3.82 0 0 1 5.218 1.398zm-17.586 8.34a1.75 1.75 0 0 1-.64 2.391 1.75 1.75 0 0 1-2.39-.64 1.75 1.75 0 0 1 .64-2.39 1.75 1.75 0 0 1 2.39.64zm6.904 2.158a4.01 4.01 0 0 1-1.468 5.478A4.01 4.01 0 0 1 265.5 93.6a4.01 4.01 0 0 1 1.468-5.478 4.01 4.01 0 0 1 5.478 1.468zm-7.969-1.542a.52.52 0 0 1-.19.71.52.52 0 0 1-.71-.19.52.52 0 0 1 .19-.71.52.52 0 0 1 .71.19zm18.365-8.79a3.49 3.49 0 0 1-1.277 4.767 3.49 3.49 0 0 1-4.768-1.278 3.49 3.49 0 0 1 1.278-4.767 3.49 3.49 0 0 1 4.767 1.277zM270.281 90.84a1.51 1.51 0 0 1-.553 2.063 1.51 1.51 0 0 1-2.063-.553 1.51 1.51 0 0 1 .553-2.063 1.51 1.51 0 0 1 2.063.553zm43.625 170.962a2.08 2.08 0 0 1-.761 2.84 2.08 2.08 0 0 1-2.841-.76 2.08 2.08 0 0 1 .76-2.842 2.08 2.08 0 0 1 2.842.762zm-24.614-47.473a1.37 1.37 0 0 1-.501 1.871 1.37 1.37 0 0 1-1.872-.501 1.37 1.37 0 0 1 .502-1.872 1.37 1.37 0 0 1 1.871.502zm45.46 45.818a3.4 3.4 0 0 1-1.245 4.645 3.4 3.4 0 0 1-4.644-1.245 3.4 3.4 0 0 1 1.244-4.644 3.4 3.4 0 0 1 4.645 1.244zm-19.945 1.135a3.12 3.12 0 0 1-1.142 4.262 3.12 3.12 0 0 1-4.262-1.142 3.12 3.12 0 0 1 1.142-4.262 3.12 3.12 0 0 1 4.262 1.142zm-47.375-52.577a3.49 3.49 0 0 1-1.278 4.768 3.49 3.49 0 0 1-4.767-1.278 3.49 3.49 0 0 1 1.277-4.767 3.49 3.49 0 0 1 4.768 1.277zm-28.067-66.033a3.68 3.68 0 0 1-1.347 5.027 3.68 3.68 0 0 1-5.027-1.347 3.68 3.68 0 0 1 1.347-5.027 3.68 3.68 0 0 1 5.027 1.347zm60.018 41.894a2.46 2.46 0 0 1-.9 3.36 2.46 2.46 0 0 1-3.361-.9 2.46 2.46 0 0 1 .9-3.36 2.46 2.46 0 0 1 3.36.9zm-34.445 25.58a.61.61 0 0 1-.224.833.61.61 0 0 1-.833-.224.61.61 0 0 1 .223-.833.61.61 0 0 1 .834.223zm1.151-.666a1.94 1.94 45 0 1-.71 2.65 1.94 1.94 45 0 1-2.65-.71 1.94 1.94 45 0 1 .71-2.65 1.94 1.94 45 0 1 2.65.71zm-22.146-.678a4.01 4.01 0 0 1-1.467 5.478 4.01 4.01 0 0 1-5.478-1.468 4.01 4.01 0 0 1 1.468-5.477 4.01 4.01 0 0 1 5.477 1.467zm-4.206-66.345a4.11 4.11 0 0 1-1.504 5.614 4.11 4.11 0 0 1-5.615-1.504 4.11 4.11 0 0 1 1.505-5.614 4.11 4.11 0 0 1 5.614 1.504zm27.539 66.338a3.31 3.31 0 0 1-1.212 4.522 3.31 3.31 0 0 1-4.521-1.212 3.31 3.31 0 0 1 1.211-4.521 3.31 3.31 0 0 1 4.522 1.211zm21.644 5.749a.94.94 45 0 1-.344 1.284.94.94 45 0 1-1.284-.344.94.94 45 0 1 .344-1.285.94.94 45 0 1 1.284.345zm-21.038-6.099a4.01 4.01 0 0 1-1.468 5.478 4.01 4.01 0 0 1-5.477-1.468 4.01 4.01 0 0 1 1.467-5.477 4.01 4.01 0 0 1 5.478 1.467zm30.184-23.12a.94.94 45 0 1-.344 1.285.94.94 45 0 1-1.284-.344.94.94 45 0 1 .344-1.284.94.94 45 0 1 1.284.344zm-6.566 27.729a3.92 3.92 45 0 1-1.434 5.354 3.92 3.92 45 0 1-5.355-1.434 3.92 3.92 45 0 1 1.435-5.355 3.92 3.92 45 0 1 5.354 1.435zm-40.407 52.693a2.08 2.08 0 0 1-.76 2.84 2.08 2.08 0 0 1-2.842-.76 2.08 2.08 0 0 1 .761-2.842 2.08 2.08 0 0 1 2.841.762zm4.71-52.522a1.98 1.98 0 0 1-.724 2.705 1.98 1.98 0 0 1-2.705-.725 1.98 1.98 0 0 1 .725-2.705 1.98 1.98 0 0 1 2.705.725zm10.408-3.815a2.08 2.08 0 0 1-.762 2.842 2.08 2.08 0 0 1-2.841-.762 2.08 2.08 0 0 1 .761-2.841 2.08 2.08 0 0 1 2.842.761zm-14.667 56.077a2.6 2.6 0 0 1-.952 3.551 2.6 2.6 0 0 1-3.552-.951 2.6 2.6 0 0 1 .952-3.552 2.6 2.6 0 0 1 3.552.952zm-1.966 1.135a.33.33 0 0 1-.12.45.33.33 0 0 1-.452-.12.33.33 0 0 1 .121-.451.33.33 0 0 1 .45.12zm15.732-56.692a1.04 1.04 0 0 1-.38 1.421 1.04 1.04 0 0 1-1.421-.38 1.04 1.04 0 0 1 .38-1.421 1.04 1.04 0 0 1 1.421.38zm23.774 4.519a1.13 1.13 0 0 1-.413 1.543 1.13 1.13 0 0 1-1.544-.413 1.13 1.13 0 0 1 .414-1.544 1.13 1.13 0 0 1 1.543.414zm-37.133 50.803a3.07 3.07 45 0 1-1.124 4.193 3.07 3.07 45 0 1-4.194-1.123 3.07 3.07 45 0 1 1.124-4.194 3.07 3.07 45 0 1 4.194 1.124zm2.346-51.157a.24.24 0 0 1-.088.328.24.24 0 0 1-.328-.088.24.24 0 0 1 .088-.328.24.24 0 0 1 .328.088zm-12.597-3.998a1.42 1.42 0 0 1-.52 1.94 1.42 1.42 0 0 1-1.939-.52 1.42 1.42 0 0 1 .52-1.94 1.42 1.42 0 0 1 1.94.52zm26.676-1.937a4.58 4.58 0 0 1-1.677 6.257 4.58 4.58 0 0 1-6.256-1.677 4.58 4.58 0 0 1 1.676-6.256 4.58 4.58 0 0 1 6.257 1.676zm-12.572 5.065a1.98 1.98 0 0 1-.725 2.705 1.98 1.98 0 0 1-2.705-.725 1.98 1.98 0 0 1 .725-2.705 1.98 1.98 0 0 1 2.705.725zm-22.652 67.707a3.45 3.45 0 0 1-1.263 4.712 3.45 3.45 0 0 1-4.712-1.262 3.45 3.45 0 0 1 1.262-4.713 3.45 3.45 0 0 1 4.713 1.263zm16.876-14.57a.85.85 0 0 1-.31 1.16.85.85 0 0 1-1.162-.31.85.85 0 0 1 .311-1.162.85.85 0 0 1 1.161.312zm17.82 8.464a.8.8 0 0 1-.293 1.093.8.8 0 0 1-1.093-.293.8.8 0 0 1 .293-1.093.8.8 0 0 1 1.093.293zm-34.08 5.75a4.16 4.16 0 0 1-1.524 5.683 4.16 4.16 0 0 1-5.682-1.522 4.16 4.16 0 0 1 1.522-5.683 4.16 4.16 0 0 1 5.683 1.523zm32.302 37.33a3.49 3.49 0 0 1-1.278 4.768 3.49 3.49 0 0 1-4.767-1.278 3.49 3.49 0 0 1 1.277-4.767 3.49 3.49 0 0 1 4.768 1.277zm-34.433-36.1a1.7 1.7 0 0 1-.623 2.323 1.7 1.7 0 0 1-2.322-.622 1.7 1.7 0 0 1 .622-2.323 1.7 1.7 0 0 1 2.323.623zm38.42-8.255a3.35 3.35 0 0 1-1.227 4.576 3.35 3.35 0 0 1-4.576-1.226 3.35 3.35 0 0 1 1.226-4.576 3.35 3.35 0 0 1 4.576 1.226zm-3.658 44.165a3.87 3.87 0 0 1-1.417 5.287 3.87 3.87 0 0 1-5.286-1.417 3.87 3.87 0 0 1 1.416-5.286 3.87 3.87 0 0 1 5.287 1.416zm-35.048-35.744a1.37 1.37 0 0 1-.502 1.871 1.37 1.37 0 0 1-1.871-.501 1.37 1.37 0 0 1 .501-1.872 1.37 1.37 0 0 1 1.872.502zm-14.183-3.486a3.12 3.12 0 0 1-1.142 4.262 3.12 3.12 0 0 1-4.262-1.142 3.12 3.12 0 0 1 1.142-4.262 3.12 3.12 0 0 1 4.262 1.142zm34.662-13.164a2.93 2.93 0 0 1-1.073 4.002 2.93 2.93 0 0 1-4.002-1.072 2.93 2.93 0 0 1 1.072-4.003 2.93 2.93 0 0 1 4.003 1.073zm-20.436 16.625a1.42 1.42 0 0 1-.52 1.94 1.42 1.42 0 0 1-1.94-.52 1.42 1.42 0 0 1 .52-1.94 1.42 1.42 0 0 1 1.94.52zm-31.178 35.437a1.65 1.65 0 0 1-.604 2.253 1.65 1.65 0 0 1-2.254-.603 1.65 1.65 0 0 1 .604-2.254 1.65 1.65 0 0 1 2.254.604zm30.52-35.057a.66.66 0 0 1-.242.901.66.66 0 0 1-.901-.241.66.66 0 0 1 .241-.902.66.66 0 0 1 .902.242zm34.563 36.024a2.6 2.6 0 0 1-.952 3.552 2.6 2.6 0 0 1-3.551-.952 2.6 2.6 0 0 1 .951-3.551 2.6 2.6 0 0 1 3.552.951zm-63.49-1.887a3.49 3.49 0 0 1-1.277 4.767 3.49 3.49 0 0 1-4.768-1.277 3.49 3.49 0 0 1 1.278-4.768 3.49 3.49 0 0 1 4.767 1.278zm-1.264.73a2.03 2.03 0 0 1-.743 2.773 2.03 2.03 0 0 1-2.773-.743 2.03 2.03 0 0 1 .743-2.773 2.03 2.03 0 0 1 2.773.743zm15.437-38.023a1.75 1.75 0 0 1-.641 2.39 1.75 1.75 0 0 1-2.39-.64 1.75 1.75 0 0 1 .64-2.39 1.75 1.75 0 0 1 2.39.64zm16.469 2.166a2.64 2.64 0 0 1-.967 3.606 2.64 2.64 0 0 1-3.606-.966 2.64 2.64 0 0 1 .966-3.607 2.64 2.64 0 0 1 3.607.967zm-29.62 34.537a4.67 4.67 0 0 1-1.71 6.379 4.67 4.67 0 0 1-6.379-1.71 4.67 4.67 0 0 1 1.71-6.379 4.67 4.67 0 0 1 6.379 1.71zm192.812-111.32a3.97 3.97 0 0 1-1.453 5.423 3.97 3.97 0 0 1-5.423-1.453 3.97 3.97 0 0 1 1.453-5.424 3.97 3.97 0 0 1 5.423 1.454zm-10.211-23.527a3.97 3.97 0 0 1-1.453 5.423 3.97 3.97 0 0 1-5.424-1.453 3.97 3.97 0 0 1 1.454-5.423 3.97 3.97 0 0 1 5.423 1.453zm-9.673-40.373a1.04 1.04 0 0 1-.38 1.421 1.04 1.04 0 0 1-1.42-.38 1.04 1.04 0 0 1 .38-1.421 1.04 1.04 0 0 1 1.42.38zm18.983 64.42a2.93 2.93 0 0 1-1.072 4.002 2.93 2.93 0 0 1-4.003-1.072A2.93 2.93 0 0 1 390.74 204a2.93 2.93 0 0 1 4.002 1.073zm-9.561-23.902a3.68 3.68 0 0 1-1.347 5.027 3.68 3.68 0 0 1-5.027-1.347 3.68 3.68 0 0 1 1.347-5.027 3.68 3.68 0 0 1 5.027 1.347zm-85.227 3.064a3.12 3.12 0 0 1-1.142 4.262 3.12 3.12 0 0 1-4.262-1.142 3.12 3.12 0 0 1 1.142-4.262 3.12 3.12 0 0 1 4.262 1.142zm76.212-43.817a1.51 1.51 0 0 1-.552 2.063 1.51 1.51 0 0 1-2.063-.553 1.51 1.51 0 0 1 .553-2.062 1.51 1.51 0 0 1 2.062.552zm9.387 40.538a4.11 4.11 0 0 1-1.504 5.615 4.11 4.11 0 0 1-5.615-1.505 4.11 4.11 0 0 1 1.505-5.614 4.11 4.11 0 0 1 5.614 1.504zm-7.22 87.156a3.59 3.59 0 0 1-1.313 4.905 3.59 3.59 0 0 1-4.904-1.315 3.59 3.59 0 0 1 1.314-4.904 3.59 3.59 0 0 1 4.904 1.314zm7.738-31.498a3.4 3.4 0 0 1-1.244 4.644 3.4 3.4 0 0 1-4.645-1.244 3.4 3.4 0 0 1 1.245-4.645 3.4 3.4 0 0 1 4.644 1.245zm7.77-31.021a1.89 1.89 0 0 1-.691 2.581 1.89 1.89 0 0 1-2.582-.691 1.89 1.89 0 0 1 .692-2.582 1.89 1.89 0 0 1 2.582.692zm-16.122 62.874a2.88 2.88 45 0 1-1.054 3.935 2.88 2.88 45 0 1-3.934-1.055 2.88 2.88 45 0 1 1.054-3.934 2.88 2.88 45 0 1 3.934 1.054zm-111.553 64.405a4.63 4.63 0 0 1-1.694 6.325 4.63 4.63 0 0 1-6.325-1.695 4.63 4.63 0 0 1 1.695-6.324 4.63 4.63 0 0 1 6.324 1.694zm8.9-11.026a4.53 4.53 0 0 1-1.659 6.188 4.53 4.53 0 0 1-6.188-1.658 4.53 4.53 0 0 1 1.658-6.188 4.53 4.53 0 0 1 6.188 1.658zm42.173-8.853a3.31 3.31 0 0 1-1.212 4.521 3.31 3.31 0 0 1-4.521-1.211 3.31 3.31 0 0 1 1.211-4.522 3.31 3.31 0 0 1 4.522 1.212zm-53.324 21.18a2.03 2.03 0 0 1-.743 2.773 2.03 2.03 0 0 1-2.773-.744 2.03 2.03 0 0 1 .743-2.773 2.03 2.03 0 0 1 2.773.743zm7.885-10.442a.76.76 0 0 1-.278 1.038.76.76 0 0 1-1.038-.278.76.76 0 0 1 .278-1.038.76.76 0 0 1 1.038.278zm-5.15-6.16a4.16 4.16 0 0 1-1.523 5.682 4.16 4.16 0 0 1-5.682-1.523 4.16 4.16 0 0 1 1.522-5.682 4.16 4.16 0 0 1 5.683 1.522zm42.292-18.61a4.16 4.16 0 0 1-1.523 5.683 4.16 4.16 0 0 1-5.683-1.522 4.16 4.16 0 0 1 1.523-5.683 4.16 4.16 0 0 1 5.683 1.523zM273.8 322.578a3.07 3.07 45 0 1-1.124 4.193 3.07 3.07 45 0 1-4.194-1.123 3.07 3.07 45 0 1 1.124-4.194 3.07 3.07 45 0 1 4.194 1.124zm.571-.33a3.73 3.73 45 0 1-1.365 5.095 3.73 3.73 45 0 1-5.095-1.365 3.73 3.73 45 0 1 1.365-5.096 3.73 3.73 45 0 1 5.095 1.366zm34.12-23.024a3.64 3.64 0 0 1-1.333 4.972 3.64 3.64 0 0 1-4.972-1.332 3.64 3.64 0 0 1 1.332-4.973 3.64 3.64 0 0 1 4.972 1.333zm7.067 14.74a1.37 1.37 0 0 1-.502 1.872 1.37 1.37 0 0 1-1.871-.501 1.37 1.37 0 0 1 .501-1.872 1.37 1.37 0 0 1 1.872.502zm-40.45 7.859a4.58 4.58 0 0 1-1.677 6.256 4.58 4.58 0 0 1-6.256-1.676 4.58 4.58 0 0 1 1.676-6.257 4.58 4.58 0 0 1 6.256 1.677zm-10.25 11.806a3.12 3.12 0 0 1-1.142 4.262 3.12 3.12 0 0 1-4.262-1.142 3.12 3.12 0 0 1 1.142-4.262 3.12 3.12 0 0 1 4.262 1.142zm.847-15.512a3.07 3.07 45 0 1-1.124 4.194 3.07 3.07 45 0 1-4.193-1.124 3.07 3.07 45 0 1 1.123-4.193 3.07 3.07 45 0 1 4.194 1.123zm8.424 4.27a3.45 3.45 0 0 1-1.263 4.714 3.45 3.45 0 0 1-4.713-1.263 3.45 3.45 0 0 1 1.263-4.713 3.45 3.45 0 0 1 4.713 1.263zm-9.523 11.387a2.83 2.83 0 0 1-1.035 3.866 2.83 2.83 0 0 1-3.866-1.036 2.83 2.83 0 0 1 1.036-3.865 2.83 2.83 0 0 1 3.865 1.035zm74.987-34.899a.94.94 45 0 1-.344 1.285.94.94 45 0 1-1.284-.345.94.94 45 0 1 .344-1.284.94.94 45 0 1 1.284.344zm38.575-30.666a3.4 3.4 0 0 1-1.244 4.645 3.4 3.4 0 0 1-4.645-1.245 3.4 3.4 0 0 1 1.245-4.644 3.4 3.4 0 0 1 4.644 1.244zm-38.081 30.381a1.51 1.51 0 0 1-.553 2.063 1.51 1.51 0 0 1-2.063-.553 1.51 1.51 0 0 1 .553-2.062 1.51 1.51 0 0 1 2.063.552zm-4.392-38.987a4.49 4.49 0 0 1-1.644 6.134 4.49 4.49 0 0 1-6.133-1.644 4.49 4.49 0 0 1 1.643-6.133 4.49 4.49 0 0 1 6.134 1.643zm43.088 8.251a4.11 4.11 0 0 1-1.504 5.615 4.11 4.11 0 0 1-5.615-1.505 4.11 4.11 0 0 1 1.505-5.614 4.11 4.11 0 0 1 5.614 1.504zM341.11 298a2.69 2.69 45 0 1-.985 3.675 2.69 2.69 45 0 1-3.674-.985 2.69 2.69 45 0 1 .984-3.674 2.69 2.69 45 0 1 3.675.984zM237.278 143.877a1.27 1.27 0 0 1-.465 1.735 1.27 1.27 0 0 1-1.735-.465 1.27 1.27 0 0 1 .465-1.735 1.27 1.27 0 0 1 1.735.465zm-34.801-4.537a1.65 1.65 0 0 1-.604 2.254 1.65 1.65 0 0 1-2.254-.604 1.65 1.65 0 0 1 .604-2.254 1.65 1.65 0 0 1 2.254.604zm50.03-29.185a2.22 2.22 0 0 1-.813 3.032 2.22 2.22 0 0 1-3.032-.812 2.22 2.22 0 0 1 .812-3.033 2.22 2.22 0 0 1 3.033.813zm-13.056 32.467a3.78 3.78 0 0 1-1.383 5.164 3.78 3.78 0 0 1-5.164-1.384 3.78 3.78 0 0 1 1.384-5.164 3.78 3.78 0 0 1 5.163 1.384zm.572-.33a4.44 4.44 0 0 1-1.625 6.065 4.44 4.44 0 0 1-6.065-1.625 4.44 4.44 0 0 1 1.625-6.065 4.44 4.44 0 0 1 6.065 1.625zm13.791-32.892a3.73 3.73 45 0 1-1.365 5.095 3.73 3.73 45 0 1-5.095-1.365 3.73 3.73 45 0 1 1.365-5.096 3.73 3.73 45 0 1 5.095 1.366zm44.581 75.736a1.32 1.32 0 0 1-.483 1.803 1.32 1.32 0 0 1-1.803-.483 1.32 1.32 0 0 1 .483-1.803 1.32 1.32 0 0 1 1.803.483zm-59.888-41.969a2.69 2.69 45 0 1-.984 3.675 2.69 2.69 45 0 1-3.675-.985 2.69 2.69 45 0 1 .985-3.675 2.69 2.69 45 0 1 3.674.985zm67.897 157.26a1.23 1.23 0 0 1-.45 1.68 1.23 1.23 0 0 1-1.68-.45 1.23 1.23 0 0 1 .45-1.68 1.23 1.23 0 0 1 1.68.45zm-35.49-27.371a4.34 4.34 0 0 1-1.59 5.928 4.34 4.34 0 0 1-5.928-1.588 4.34 4.34 0 0 1 1.589-5.929 4.34 4.34 0 0 1 5.928 1.589zm44.084 11.256a2.5 2.5 0 0 1-.915 3.415 2.5 2.5 0 0 1-3.416-.915 2.5 2.5 0 0 1 .916-3.415 2.5 2.5 0 0 1 3.415.915zm-7.78 15.646a2.17 2.17 0 0 1-.794 2.964 2.17 2.17 0 0 1-2.964-.794 2.17 2.17 0 0 1 .794-2.965 2.17 2.17 0 0 1 2.964.795zm-38.184-25.817a2.17 2.17 0 0 1-.794 2.964 2.17 2.17 0 0 1-2.964-.794 2.17 2.17 0 0 1 .794-2.964 2.17 2.17 0 0 1 2.964.794zm-18.928-7.824a.94.94 45 0 1-.344 1.284.94.94 45 0 1-1.284-.344.94.94 45 0 1 .344-1.284.94.94 45 0 1 1.284.344zm39.264-52.033a1.46 1.46 0 0 1-.534 1.994 1.46 1.46 0 0 1-1.995-.534 1.46 1.46 0 0 1 .535-1.995 1.46 1.46 0 0 1 1.994.535zm-21.6 60.587a.71.71 0 0 1-.26.97.71.71 0 0 1-.97-.26.71.71 0 0 1 .26-.97.71.71 0 0 1 .97.26zm44.863-12.334a.61.61 0 0 1-.223.833.61.61 0 0 1-.833-.223.61.61 0 0 1 .223-.834.61.61 0 0 1 .833.224zm-42.49 10.964a3.45 3.45 0 0 1-1.263 4.713 3.45 3.45 0 0 1-4.713-1.263 3.45 3.45 0 0 1 1.263-4.713 3.45 3.45 0 0 1 4.713 1.263zm21.314-60.422a3.87 3.87 0 0 1-1.416 5.286 3.87 3.87 0 0 1-5.287-1.416 3.87 3.87 0 0 1 1.417-5.287 3.87 3.87 0 0 1 5.286 1.417zm23.307 48.228a3.07 3.07 45 0 1-1.124 4.193 3.07 3.07 45 0 1-4.194-1.123 3.07 3.07 45 0 1 1.124-4.194 3.07 3.07 45 0 1 4.194 1.124zm-1.195 23.83a.85.85 0 0 1-.311 1.16.85.85 0 0 1-1.162-.31.85.85 0 0 1 .312-1.161.85.85 0 0 1 1.16.31zM267.89 274.8a.85.85 0 0 1-.311 1.161.85.85 0 0 1-1.161-.311.85.85 0 0 1 .31-1.161.85.85 0 0 1 1.162.31zm47.08-13.614a3.31 3.31 0 0 1-1.211 4.521 3.31 3.31 0 0 1-4.522-1.211 3.31 3.31 0 0 1 1.212-4.522 3.31 3.31 0 0 1 4.521 1.212zm-.787 23.595a1.56 1.56 0 0 1-.571 2.13 1.56 1.56 0 0 1-2.131-.57 1.56 1.56 0 0 1 .57-2.131 1.56 1.56 0 0 1 2.132.57zM301.01 183.626a4.34 4.34 0 0 1-1.589 5.928 4.34 4.34 0 0 1-5.928-1.588 4.34 4.34 0 0 1 1.588-5.929 4.34 4.34 0 0 1 5.929 1.589zm49.237-65.239a2.46 2.46 0 0 1-.9 3.36 2.46 2.46 0 0 1-3.361-.9 2.46 2.46 0 0 1 .9-3.36 2.46 2.46 0 0 1 3.36.9zm26.94 21.442a2.69 2.69 45 0 1-.984 3.675 2.69 2.69 45 0 1-3.675-.985 2.69 2.69 45 0 1 .985-3.674 2.69 2.69 45 0 1 3.674.984zm-79.278 45.587a.76.76 0 0 1-.278 1.038.76.76 0 0 1-1.038-.278.76.76 0 0 1 .278-1.038.76.76 0 0 1 1.038.278zm-144.615 66.439a3.21 3.21 0 0 1-1.175 4.385 3.21 3.21 0 0 1-4.385-1.175 3.21 3.21 0 0 1 1.175-4.385 3.21 3.21 0 0 1 4.385 1.175zm-18.718-47.321a1.23 1.23 0 0 1-.45 1.68 1.23 1.23 0 0 1-1.68-.45 1.23 1.23 0 0 1 .45-1.68 1.23 1.23 0 0 1 1.68.45zm50.57 21.23a2.27 2.27 0 0 1-.831 3.1 2.27 2.27 0 0 1-3.101-.83 2.27 2.27 0 0 1 .83-3.102 2.27 2.27 0 0 1 3.102.831zm-30.622 25.38a4.63 4.63 0 0 1-1.695 6.325 4.63 4.63 0 0 1-6.325-1.694 4.63 4.63 0 0 1 1.695-6.325 4.63 4.63 0 0 1 6.325 1.695zm36.333-11.589a3.26 3.26 45 0 1-1.193 4.454 3.26 3.26 45 0 1-4.454-1.194 3.26 3.26 45 0 1 1.194-4.453 3.26 3.26 45 0 1 4.453 1.193zm-38.013 12.56a2.69 2.69 45 0 1-.985 3.674 2.69 2.69 45 0 1-3.675-.984 2.69 2.69 45 0 1 .985-3.675 2.69 2.69 45 0 1 3.675.985zm33.687-27.152a3.87 3.87 0 0 1-1.416 5.287 3.87 3.87 0 0 1-5.287-1.417 3.87 3.87 0 0 1 1.417-5.286 3.87 3.87 0 0 1 5.286 1.416zm3.217 15.232a1.98 1.98 0 0 1-.724 2.705 1.98 1.98 0 0 1-2.705-.725 1.98 1.98 0 0 1 .725-2.704 1.98 1.98 0 0 1 2.704.724zm-38.403 30.024a.57.57 0 0 1-.208.779.57.57 0 0 1-.779-.209.57.57 0 0 1 .209-.778.57.57 0 0 1 .778.208zm1.499-18.104a2.69 2.69 45 0 1-.985 3.674 2.69 2.69 45 0 1-3.675-.984 2.69 2.69 45 0 1 .985-3.675 2.69 2.69 45 0 1 3.675.985zm35.519-11.12a.38.38 0 0 1-.14.52.38.38 0 0 1-.519-.14.38.38 0 0 1 .14-.519.38.38 0 0 1 .519.14zm-33.909 27.43a4.16 4.16 0 0 1-1.522 5.682 4.16 4.16 0 0 1-5.683-1.523 4.16 4.16 0 0 1 1.523-5.682 4.16 4.16 0 0 1 5.682 1.522zm-2.754-15.65a1.37 1.37 0 0 1-.501 1.871 1.37 1.37 0 0 1-1.871-.501 1.37 1.37 0 0 1 .501-1.872 1.37 1.37 0 0 1 1.871.502zm2.668 15.7a4.06 4.06 0 0 1-1.486 5.545 4.06 4.06 0 0 1-5.546-1.486 4.06 4.06 0 0 1 1.486-5.546 4.06 4.06 0 0 1 5.546 1.486zm-18.034-64.956a3.26 3.26 45 0 1-1.193 4.453 3.26 3.26 45 0 1-4.454-1.193 3.26 3.26 45 0 1 1.194-4.453 3.26 3.26 45 0 1 4.453 1.193zm17.982 47.746a4.39 4.39 0 0 1-1.607 5.996 4.39 4.39 0 0 1-5.997-1.606 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.997 1.607zm41.456 46.964a4.3 4.3 45 0 1-1.574 5.874 4.3 4.3 45 0 1-5.874-1.574 4.3 4.3 45 0 1 1.574-5.874 4.3 4.3 45 0 1 5.874 1.574zm-40.997-29.99a4.53 4.53 0 0 1-1.658 6.189 4.53 4.53 0 0 1-6.188-1.659 4.53 4.53 0 0 1 1.658-6.188 4.53 4.53 0 0 1 6.188 1.658zm35.017-28.069a2.03 2.03 0 0 1-.743 2.773 2.03 2.03 0 0 1-2.773-.743 2.03 2.03 0 0 1 .743-2.773 2.03 2.03 0 0 1 2.773.743zm5.278 58.464a3.49 3.49 0 0 1-1.277 4.767 3.49 3.49 0 0 1-4.767-1.277 3.49 3.49 0 0 1 1.277-4.768 3.49 3.49 0 0 1 4.767 1.278zm69.865 34.95a3.21 3.21 0 0 1-1.174 4.384 3.21 3.21 0 0 1-4.385-1.175 3.21 3.21 0 0 1 1.175-4.384 3.21 3.21 0 0 1 4.385 1.174zm50.294-19.43a.99.99 45 0 1-.363 1.352.99.99 45 0 1-1.352-.362.99.99 45 0 1 .362-1.353.99.99 45 0 1 1.353.363zm-51.601 20.184a1.7 1.7 0 0 1-.622 2.323 1.7 1.7 0 0 1-2.323-.623 1.7 1.7 0 0 1 .623-2.322 1.7 1.7 0 0 1 2.322.622zm.726-15.442a1.51 1.51 0 0 1-.553 2.063 1.51 1.51 0 0 1-2.062-.553 1.51 1.51 0 0 1 .552-2.062 1.51 1.51 0 0 1 2.063.552zm4.48-44.64a1.94 1.94 45 0 1-.71 2.65 1.94 1.94 45 0 1-2.65-.71 1.94 1.94 45 0 1 .71-2.65 1.94 1.94 45 0 1 2.65.71zm37.361 26.292a.99.99 45 0 1-.362 1.352.99.99 45 0 1-1.352-.362.99.99 45 0 1 .362-1.353.99.99 45 0 1 1.352.363zm-42.005 18.443a1.32 1.32 0 0 1-.484 1.804 1.32 1.32 0 0 1-1.803-.484 1.32 1.32 0 0 1 .483-1.803 1.32 1.32 0 0 1 1.804.483zm-11.805-208.766a2.08 2.08 0 0 1-.762 2.84 2.08 2.08 0 0 1-2.841-.76 2.08 2.08 0 0 1 .761-2.842 2.08 2.08 0 0 1 2.842.762zm17.159-18.96a.66.66 0 0 1-.242.901.66.66 0 0 1-.902-.241.66.66 0 0 1 .242-.902.66.66 0 0 1 .902.242zm34.676 2.461a1.09 1.09 0 0 1-.399 1.49 1.09 1.09 0 0 1-1.489-.4 1.09 1.09 0 0 1 .4-1.489 1.09 1.09 0 0 1 1.488.4zm-52.615 16.949a1.18 1.18 0 0 1-.432 1.611 1.18 1.18 0 0 1-1.612-.431 1.18 1.18 0 0 1 .432-1.612 1.18 1.18 0 0 1 1.612.432zm-109.136 73.31a1.09 1.09 0 0 1-.4 1.488 1.09 1.09 0 0 1-1.488-.399 1.09 1.09 0 0 1 .399-1.488 1.09 1.09 0 0 1 1.489.398zm-6.871 19.959a2.41 2.41 0 0 1-.882 3.292 2.41 2.41 0 0 1-3.292-.882 2.41 2.41 0 0 1 .882-3.292 2.41 2.41 0 0 1 3.292.882zm18.447-64.229a4.11 4.11 0 0 1-1.504 5.615 4.11 4.11 0 0 1-5.615-1.505 4.11 4.11 0 0 1 1.505-5.614 4.11 4.11 0 0 1 5.614 1.504zm-12.191 44.624a.38.38 0 0 1-.14.52.38.38 0 0 1-.519-.14.38.38 0 0 1 .14-.519.38.38 0 0 1 .519.14zm.199-.115a.61.61 0 0 1-.223.834.61.61 0 0 1-.834-.224.61.61 0 0 1 .224-.833.61.61 0 0 1 .833.223zm12.607-38.744a.66.66 0 0 1-.242.902.66.66 0 0 1-.901-.242.66.66 0 0 1 .241-.901.66.66 0 0 1 .902.241zm41.78 22.806a3.07 3.07 45 0 1-1.124 4.193 3.07 3.07 45 0 1-4.194-1.123 3.07 3.07 45 0 1 1.124-4.194 3.07 3.07 45 0 1 4.194 1.124zm-54.223 15.843a.8.8 0 0 1-.292 1.093.8.8 0 0 1-1.093-.293.8.8 0 0 1 .293-1.092.8.8 0 0 1 1.092.292zm46.628 36.681a3.45 3.45 0 0 1-1.263 4.713 3.45 3.45 0 0 1-4.713-1.263 3.45 3.45 0 0 1 1.263-4.713 3.45 3.45 0 0 1 4.713 1.263zm-46.749-36.61a.66.66 0 0 1-.241.901.66.66 0 0 1-.902-.242.66.66 0 0 1 .242-.901.66.66 0 0 1 .901.241zm55.123-16.364a3.97 3.97 0 0 1-1.453 5.423 3.97 3.97 0 0 1-5.423-1.453 3.97 3.97 0 0 1 1.453-5.424 3.97 3.97 0 0 1 5.423 1.454zm-10.462 54.18a1.04 1.04 0 0 1-.38 1.42 1.04 1.04 0 0 1-1.42-.38 1.04 1.04 0 0 1 .38-1.422 1.04 1.04 0 0 1 1.42.381zm-52.103-17.527a1.32 1.32 0 0 1-.483 1.803 1.32 1.32 0 0 1-1.804-.483 1.32 1.32 0 0 1 .484-1.803 1.32 1.32 0 0 1 1.803.483zm9.737-21.615a3.31 3.31 0 0 1-1.211 4.522 3.31 3.31 0 0 1-4.522-1.212 3.31 3.31 0 0 1 1.212-4.521 3.31 3.31 0 0 1 4.521 1.211zm43.146 38.691a1.94 1.94 45 0 1-.71 2.65 1.94 1.94 45 0 1-2.65-.71 1.94 1.94 45 0 1 .71-2.65 1.94 1.94 45 0 1 2.65.71zm-52.39-17.361a1.89 1.89 0 0 1-.691 2.582 1.89 1.89 0 0 1-2.582-.692 1.89 1.89 0 0 1 .692-2.582 1.89 1.89 0 0 1 2.581.692zm8.014-20.62a1.89 1.89 0 0 1-.691 2.582 1.89 1.89 0 0 1-2.582-.692 1.89 1.89 0 0 1 .692-2.581 1.89 1.89 0 0 1 2.581.691zm8.06-42.239a.85.85 0 0 1-.31 1.161.85.85 0 0 1-1.162-.31.85.85 0 0 1 .312-1.162.85.85 0 0 1 1.16.311zm3.439 4.135a.66.66 0 0 1-.242.902.66.66 0 0 1-.901-.242.66.66 0 0 1 .241-.901.66.66 0 0 1 .902.241zm-12.278 38.554a.99.99 45 0 1-.362 1.353.99.99 45 0 1-1.353-.363.99.99 45 0 1 .363-1.352.99.99 45 0 1 1.352.362zm209.58-66.637a4.44 4.44 0 0 1-1.626 6.065 4.44 4.44 0 0 1-6.065-1.625 4.44 4.44 0 0 1 1.625-6.065 4.44 4.44 0 0 1 6.065 1.625zm-23.067-25.053a2.27 2.27 0 0 1-.831 3.1 2.27 2.27 0 0 1-3.101-.83 2.27 2.27 0 0 1 .83-3.1 2.27 2.27 0 0 1 3.102.83zm48.535 47.345a2.97 2.97 0 0 1-1.087 4.058 2.97 2.97 0 0 1-4.057-1.088 2.97 2.97 0 0 1 1.087-4.057 2.97 2.97 0 0 1 4.057 1.087zm-25.754-22.127a4.11 4.11 0 0 1-1.505 5.614 4.11 4.11 0 0 1-5.614-1.504 4.11 4.11 0 0 1 1.504-5.614 4.11 4.11 0 0 1 5.615 1.504zm-47.3-23.926a1.27 1.27 0 0 1-.465 1.735 1.27 1.27 0 0 1-1.735-.465 1.27 1.27 0 0 1 .465-1.734 1.27 1.27 0 0 1 1.735.464zm25.333-1.762a3.21 3.21 0 0 1-1.175 4.385 3.21 3.21 0 0 1-4.385-1.175 3.21 3.21 0 0 1 1.175-4.385 3.21 3.21 0 0 1 4.385 1.175zm18.65 27.603a.28.28 0 0 1-.103.382.28.28 0 0 1-.382-.102.28.28 0 0 1 .102-.383.28.28 0 0 1 .383.103zm-42.424-26.74a3.07 3.07 45 0 1-1.124 4.193 3.07 3.07 45 0 1-4.194-1.124 3.07 3.07 45 0 1 1.124-4.193 3.07 3.07 45 0 1 4.194 1.123zm-5.289 91.099a3.92 3.92 45 0 1-1.435 5.355 3.92 3.92 45 0 1-5.355-1.435 3.92 3.92 45 0 1 1.435-5.355 3.92 3.92 45 0 1 5.355 1.435zm5.618-91.29a3.45 3.45 0 0 1-1.263 4.713 3.45 3.45 0 0 1-4.713-1.263 3.45 3.45 0 0 1 1.263-4.712 3.45 3.45 0 0 1 4.713 1.262zm45.654 24.876a4.39 4.39 0 0 1-1.607 5.997 4.39 4.39 0 0 1-5.997-1.607 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.997 1.607zm-53.437 67.664a1.42 1.42 0 0 1-.52 1.94 1.42 1.42 0 0 1-1.94-.52 1.42 1.42 0 0 1 .52-1.94 1.42 1.42 0 0 1 1.94.52zm-.702.405a.61.61 0 0 1-.223.833.61.61 0 0 1-.833-.223.61.61 0 0 1 .223-.833.61.61 0 0 1 .833.223zm-46.953-74.366a.28.28 0 0 1-.103.382.28.28 0 0 1-.382-.102.28.28 0 0 1 .102-.383.28.28 0 0 1 .383.103zm54.33-17.939a2.17 2.17 0 0 1-.795 2.965 2.17 2.17 0 0 1-2.964-.795 2.17 2.17 0 0 1 .794-2.964 2.17 2.17 0 0 1 2.964.794zm-5.619 91.29a2.64 2.64 0 0 1-.966 3.606 2.64 2.64 0 0 1-3.606-.966 2.64 2.64 0 0 1 .966-3.606 2.64 2.64 0 0 1 3.606.966zm77.16 84.582a1.7 1.7 0 0 1-.623 2.323 1.7 1.7 0 0 1-2.322-.623 1.7 1.7 0 0 1 .622-2.322 1.7 1.7 0 0 1 2.322.622zm-42.353-8.676a2.93 2.93 0 0 1-1.073 4.003 2.93 2.93 0 0 1-4.002-1.073 2.93 2.93 0 0 1 1.072-4.002 2.93 2.93 0 0 1 4.003 1.072zm50.748-23.202a2.27 2.27 0 0 1-.831 3.1 2.27 2.27 0 0 1-3.101-.83 2.27 2.27 0 0 1 .83-3.101 2.27 2.27 0 0 1 3.102.83zm-9.175 32.328a.8.8 0 0 1-.293 1.093.8.8 0 0 1-1.093-.293.8.8 0 0 1 .293-1.092.8.8 0 0 1 1.093.292zm-109.64-182.5a2.6 2.6 0 0 1-.951 3.551 2.6 2.6 0 0 1-3.552-.951 2.6 2.6 0 0 1 .952-3.552 2.6 2.6 0 0 1 3.551.952zm-20.629-2.35a4.39 4.39 0 0 1-1.606 5.996 4.39 4.39 0 0 1-5.997-1.607 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.996 1.607zm21.3-10.128a3.92 3.92 45 0 1-1.434 5.355 3.92 3.92 45 0 1-5.355-1.435 3.92 3.92 45 0 1 1.435-5.354 3.92 3.92 45 0 1 5.355 1.434zm-1.693 13.068a1.42 1.42 0 0 1-.52 1.94 1.42 1.42 0 0 1-1.94-.52 1.42 1.42 0 0 1 .52-1.94 1.42 1.42 0 0 1 1.94.52zm-11.52 21.847a3.64 3.64 0 0 1-1.332 4.972 3.64 3.64 0 0 1-4.972-1.332 3.64 3.64 0 0 1 1.332-4.973 3.64 3.64 0 0 1 4.973 1.333zm11.434-21.797a1.32 1.32 0 0 1-.483 1.803 1.32 1.32 0 0 1-1.803-.483 1.32 1.32 0 0 1 .483-1.804 1.32 1.32 0 0 1 1.803.484zm4.01 3.827a.24.24 0 0 1-.087.328.24.24 0 0 1-.328-.088.24.24 0 0 1 .088-.328.24.24 0 0 1 .328.088zm-14.707 17.545a4.49 4.49 0 0 1-1.644 6.133 4.49 4.49 0 0 1-6.133-1.643 4.49 4.49 0 0 1 1.643-6.134 4.49 4.49 0 0 1 6.134 1.644zm-3.196 1.845a.8.8 0 0 1-.293 1.092.8.8 0 0 1-1.093-.292.8.8 0 0 1 .293-1.093.8.8 0 0 1 1.093.293zm-6.65-25.618a3.21 3.21 0 0 1-1.174 4.385 3.21 3.21 0 0 1-4.385-1.175 3.21 3.21 0 0 1 1.175-4.385 3.21 3.21 0 0 1 4.384 1.175zm20.993 2.14a1.84 1.84 0 0 1-.673 2.514 1.84 1.84 0 0 1-2.514-.673 1.84 1.84 0 0 1 .674-2.514 1.84 1.84 0 0 1 2.513.674zm-11.849 22.038a3.68 3.68 0 0 1-1.347 5.026 3.68 3.68 0 0 1-5.027-1.346 3.68 3.68 0 0 1 1.347-5.027 3.68 3.68 0 0 1 5.027 1.347zm-11.378 100.272a2.22 2.22 0 0 1-.812 3.033 2.22 2.22 0 0 1-3.033-.813 2.22 2.22 0 0 1 .813-3.032 2.22 2.22 0 0 1 3.032.812zm-45.216-41.836a3.92 3.92 45 0 1-1.435 5.354 3.92 3.92 45 0 1-5.355-1.434 3.92 3.92 45 0 1 1.435-5.355 3.92 3.92 45 0 1 5.355 1.435zm39.78-23.799a.9.9 0 0 1-.33 1.23.9.9 0 0 1-1.229-.33.9.9 0 0 1 .33-1.23.9.9 0 0 1 1.23.33zm5.228 65.755a1.98 1.98 0 0 1-.724 2.705 1.98 1.98 0 0 1-2.705-.725 1.98 1.98 0 0 1 .725-2.704 1.98 1.98 0 0 1 2.704.724zm-46.316-41.201a2.41 2.41 0 0 1-.882 3.292 2.41 2.41 0 0 1-3.292-.882 2.41 2.41 0 0 1 .882-3.293 2.41 2.41 0 0 1 3.292.883zm5.915-28.876a.85.85 0 0 1-.311 1.16.85.85 0 0 1-1.161-.31.85.85 0 0 1 .31-1.161.85.85 0 0 1 1.162.31zm34.844 4.512a.52.52 0 0 1-.19.71.52.52 0 0 1-.71-.19.52.52 0 0 1 .19-.71.52.52 0 0 1 .71.19zm-40.023 23.939a3.26 3.26 45 0 1-1.193 4.453 3.26 3.26 45 0 1-4.453-1.193 3.26 3.26 45 0 1 1.193-4.454 3.26 3.26 45 0 1 4.453 1.194zm-7.31 52.36a3.97 3.97 0 0 1-1.452 5.422 3.97 3.97 0 0 1-5.423-1.453 3.97 3.97 0 0 1 1.453-5.423 3.97 3.97 0 0 1 5.423 1.453zm7.596-52.525a3.59 3.59 0 0 1-1.314 4.904 3.59 3.59 0 0 1-4.904-1.314 3.59 3.59 0 0 1 1.314-4.904 3.59 3.59 0 0 1 4.904 1.314zm45.416 41.721a2.12 2.12 45 0 1-.776 2.896 2.12 2.12 45 0 1-2.896-.776 2.12 2.12 45 0 1 .776-2.895 2.12 2.12 45 0 1 2.896.775zm-54.899 11.893a1.79 1.79 0 0 1-.655 2.445 1.79 1.79 0 0 1-2.445-.655 1.79 1.79 0 0 1 .655-2.445 1.79 1.79 0 0 1 2.445.655zm63.763 44.062a2.17 2.17 0 0 1-.794 2.964 2.17 2.17 0 0 1-2.964-.794 2.17 2.17 0 0 1 .794-2.965 2.17 2.17 0 0 1 2.964.795zm-10.337-55.105a.42.42 0 0 1-.153.574.42.42 0 0 1-.574-.154.42.42 0 0 1 .154-.573.42.42 0 0 1 .573.153zm17.221 1.328a4.58 4.58 0 0 1-1.676 6.256 4.58 4.58 0 0 1-6.256-1.676 4.58 4.58 0 0 1 1.676-6.256 4.58 4.58 0 0 1 6.256 1.676zm-6.962 53.822a2.08 2.08 0 0 1-.76 2.84 2.08 2.08 0 0 1-2.842-.76 2.08 2.08 0 0 1 .761-2.842 2.08 2.08 0 0 1 2.841.762zm.936-.54a3.16 3.16 0 0 1-1.157 4.316 3.16 3.16 0 0 1-4.317-1.156 3.16 3.16 0 0 1 1.157-4.317 3.16 3.16 0 0 1 4.317 1.157zm-60.764-25.887a3.73 3.73 45 0 1-1.365 5.096 3.73 3.73 45 0 1-5.096-1.366 3.73 3.73 45 0 1 1.366-5.095 3.73 3.73 45 0 1 5.095 1.365zm50.8-29.433a1.84 1.84 0 0 1-.674 2.514 1.84 1.84 0 0 1-2.514-.674 1.84 1.84 0 0 1 .674-2.513 1.84 1.84 0 0 1 2.513.673zm7.677 56.64a.52.52 0 0 1-.19.71.52.52 0 0 1-.71-.19.52.52 0 0 1 .19-.71.52.52 0 0 1 .71.19zm-55.345-97.061a.71.71 0 0 1-.26.97.71.71 0 0 1-.97-.26.71.71 0 0 1 .26-.97.71.71 0 0 1 .97.26zm-24.543-32.55a4.53 4.53 0 0 1-1.658 6.188 4.53 4.53 0 0 1-6.188-1.658 4.53 4.53 0 0 1 1.658-6.188 4.53 4.53 0 0 1 6.188 1.658zm32.138 2.704a1.09 1.09 0 0 1-.4 1.489 1.09 1.09 0 0 1-1.488-.4 1.09 1.09 0 0 1 .399-1.488 1.09 1.09 0 0 1 1.489.399zm-6.66 29.306a1.79 1.79 0 0 1-.655 2.445 1.79 1.79 0 0 1-2.445-.655 1.79 1.79 0 0 1 .655-2.446 1.79 1.79 0 0 1 2.445.656zm-3.288 69.944a4.63 4.63 0 0 1-1.694 6.325 4.63 4.63 0 0 1-6.325-1.695 4.63 4.63 0 0 1 1.695-6.324 4.63 4.63 0 0 1 6.324 1.694zm-5.493-16.735a.8.8 0 0 1-.293 1.093.8.8 0 0 1-1.093-.293.8.8 0 0 1 .293-1.093.8.8 0 0 1 1.093.293zm56.786-12.983a3.31 3.31 0 0 1-1.211 4.522 3.31 3.31 0 0 1-4.522-1.212 3.31 3.31 0 0 1 1.212-4.521 3.31 3.31 0 0 1 4.521 1.211zM189.43 240.38a1.61 1.61 0 0 1-.59 2.2 1.61 1.61 0 0 1-2.199-.59 1.61 1.61 0 0 1 .59-2.199 1.61 1.61 0 0 1 2.199.59zm25.609 39.336a.66.66 0 0 1-.242.902.66.66 0 0 1-.902-.242.66.66 0 0 1 .242-.902.66.66 0 0 1 .902.242zm-24.838-39.78a2.5 2.5 0 0 1-.915 3.414 2.5 2.5 0 0 1-3.415-.915 2.5 2.5 0 0 1 .915-3.415 2.5 2.5 0 0 1 3.415.915zm59.5 26.616a.47.47 45 0 1-.173.642.47.47 45 0 1-.642-.172.47.47 45 0 1 .172-.642.47.47 45 0 1 .642.172zm-31.553 11.369a4.25 4.25 0 0 1-1.556 5.806 4.25 4.25 0 0 1-5.806-1.556 4.25 4.25 0 0 1 1.556-5.806 4.25 4.25 0 0 1 5.806 1.556zm-62.092-133.246a2.27 2.27 0 0 1-.83 3.101 2.27 2.27 0 0 1-3.102-.83 2.27 2.27 0 0 1 .831-3.102 2.27 2.27 0 0 1 3.101.831zm12.37-6.934a2.88 2.88 45 0 1-1.054 3.934 2.88 2.88 45 0 1-3.934-1.054 2.88 2.88 45 0 1 1.054-3.934 2.88 2.88 45 0 1 3.934 1.054zm27.072 31.09a1.98 1.98 0 0 1-.725 2.704 1.98 1.98 0 0 1-2.705-.724 1.98 1.98 0 0 1 .725-2.705 1.98 1.98 0 0 1 2.705.725zm-39.529-24.106a2.17 2.17 0 0 1-.794 2.965 2.17 2.17 0 0 1-2.964-.795 2.17 2.17 0 0 1 .794-2.964 2.17 2.17 0 0 1 2.964.794zm.494-.285a2.74 2.74 0 0 1-1.003 3.743 2.74 2.74 0 0 1-3.743-1.003 2.74 2.74 0 0 1 1.003-3.743 2.74 2.74 0 0 1 3.743 1.003zm-4.01-3.805a2.27 2.27 0 0 1-.83 3.101 2.27 2.27 0 0 1-3.101-.83 2.27 2.27 0 0 1 .83-3.102 2.27 2.27 0 0 1 3.101.831zm14.744-2.184a1.46 1.46 0 0 1-.535 1.995 1.46 1.46 0 0 1-1.994-.535 1.46 1.46 0 0 1 .534-1.994 1.46 1.46 0 0 1 1.995.534zm-10.284 5.73a3.26 3.26 45 0 1-1.193 4.453 3.26 3.26 45 0 1-4.453-1.194 3.26 3.26 45 0 1 1.193-4.453 3.26 3.26 45 0 1 4.453 1.193zm106.596 190.227a1.56 1.56 0 0 1-.571 2.131 1.56 1.56 0 0 1-2.131-.57 1.56 1.56 0 0 1 .57-2.132 1.56 1.56 0 0 1 2.132.571zm-61.898-17.83a3.26 3.26 45 0 1-1.193 4.454 3.26 3.26 45 0 1-4.454-1.193 3.26 3.26 45 0 1 1.194-4.454 3.26 3.26 45 0 1 4.453 1.194zm65.196.903a4.34 4.34 0 0 1-1.589 5.929 4.34 4.34 0 0 1-5.928-1.589 4.34 4.34 0 0 1 1.588-5.928 4.34 4.34 0 0 1 5.929 1.588zm-2.233 16.312a2.79 2.79 0 0 1-1.021 3.812 2.79 2.79 0 0 1-3.812-1.022 2.79 2.79 0 0 1 1.022-3.81 2.79 2.79 0 0 1 3.81 1.02zm-79.756-107.84a1.89 1.89 0 0 1-.691 2.582 1.89 1.89 0 0 1-2.582-.692 1.89 1.89 0 0 1 .692-2.582 1.89 1.89 0 0 1 2.581.692zm-47.833-22.81a4.01 4.01 0 0 1-1.467 5.479 4.01 4.01 0 0 1-5.478-1.468 4.01 4.01 0 0 1 1.468-5.478 4.01 4.01 0 0 1 5.477 1.468zm49.203 19.202a.38.38 0 0 1-.14.52.38.38 0 0 1-.518-.14.38.38 0 0 1 .139-.519.38.38 0 0 1 .519.14zm-.512 3.113a2.88 2.88 45 0 1-1.054 3.934 2.88 2.88 45 0 1-3.934-1.054 2.88 2.88 45 0 1 1.054-3.934 2.88 2.88 45 0 1 3.934 1.054zm4.568 14.452a2.55 2.55 0 0 1-.933 3.484 2.55 2.55 0 0 1-3.484-.934 2.55 2.55 0 0 1 .934-3.483 2.55 2.55 0 0 1 3.483.933zm-3.997-14.782a3.54 3.54 0 0 1-1.295 4.836 3.54 3.54 0 0 1-4.836-1.296 3.54 3.54 0 0 1 1.296-4.835 3.54 3.54 0 0 1 4.835 1.295zm.677-3.208a1.23 1.23 0 0 1-.45 1.68 1.23 1.23 0 0 1-1.68-.45 1.23 1.23 0 0 1 .45-1.68 1.23 1.23 0 0 1 1.68.45zm4.056 17.565a3.4 3.4 0 0 1-1.244 4.645 3.4 3.4 0 0 1-4.645-1.245 3.4 3.4 0 0 1 1.245-4.644 3.4 3.4 0 0 1 4.644 1.244zm2.049 60.329a1.13 1.13 0 0 1-.414 1.543 1.13 1.13 0 0 1-1.544-.413 1.13 1.13 0 0 1 .414-1.544 1.13 1.13 0 0 1 1.544.414zm-4.379-58.984a.71.71 0 0 1-.26.97.71.71 0 0 1-.97-.26.71.71 0 0 1 .26-.97.71.71 0 0 1 .97.26zm26.389 38.886a.66.66 0 0 1-.242.902.66.66 0 0 1-.902-.242.66.66 0 0 1 .242-.902.66.66 0 0 1 .902.242zm-18.98 18.348a4.63 4.63 0 0 1-1.694 6.324 4.63 4.63 0 0 1-6.325-1.694 4.63 4.63 0 0 1 1.695-6.325 4.63 4.63 0 0 1 6.325 1.695zm2.971 20.005a.28.28 0 0 1-.102.382.28.28 0 0 1-.383-.102.28.28 0 0 1 .103-.383.28.28 0 0 1 .382.103zm-6.452-17.995a.61.61 0 0 1-.223.833.61.61 0 0 1-.833-.223.61.61 0 0 1 .223-.834.61.61 0 0 1 .833.224zm25.69-22.223a4.39 4.39 0 0 1-1.606 5.997 4.39 4.39 0 0 1-5.997-1.607 4.39 4.39 0 0 1 1.607-5.997 4.39 4.39 0 0 1 5.997 1.607zm-16.129 38.423a3.87 3.87 0 0 1-1.416 5.286 3.87 3.87 0 0 1-5.287-1.416 3.87 3.87 0 0 1 1.417-5.287 3.87 3.87 0 0 1 5.286 1.417zm-.528.305a3.26 3.26 45 0 1-1.193 4.453 3.26 3.26 45 0 1-4.454-1.193 3.26 3.26 45 0 1 1.194-4.454 3.26 3.26 45 0 1 4.453 1.194zm-49.735-46.665a1.18 1.18 0 0 1-.432 1.612 1.18 1.18 0 0 1-1.612-.432 1.18 1.18 0 0 1 .432-1.611 1.18 1.18 0 0 1 1.612.431zm44.184 28.15a4.63 4.63 0 0 1-1.695 6.324 4.63 4.63 0 0 1-6.325-1.694 4.63 4.63 0 0 1 1.695-6.325 4.63 4.63 0 0 1 6.325 1.695zm5.837 18.35a3.59 3.59 0 0 1-1.314 4.904 3.59 3.59 0 0 1-4.904-1.314 3.59 3.59 0 0 1 1.314-4.904 3.59 3.59 0 0 1 4.904 1.314zm113.433-13.368-2.497-17.484 25.948 13.782-23.45 3.702m12.19 14.796-12.19-14.796 23.45-3.702-11.26 18.498m-13.148-3.193.957-11.603 12.192 14.796-13.149-3.193m.957-11.603-9.99-2.003 7.494-15.481 2.496 17.484m-.957 11.603-9.033-13.606 9.99 2.003-.957 11.603M241.848 86.852l-43.657 6.984 65.364-17.346-21.707 10.362m8.736 24.413-52.393-17.429 43.657-6.984 8.736 24.413m-49.536 28.9-2.857-46.329 52.393 17.429-49.536 28.9m-35.117-.984-15.444 2.59 47.704-47.935-32.26 45.345m35.117.984-35.117-.984 32.26-45.345 2.857 46.329m163.83 76.643-67.626-31.012 84.742-2.784-17.116 33.796m27.327-10.27-27.327 10.27 17.116-33.796 10.21 23.527m-9.077 31.776-18.249-21.507 27.327-10.27-9.078 31.777m-95.021-23.301 9.146-29.218 67.626 31.012-76.772-1.794m43.701 46.833 33.071-45.04 18.249 21.508-51.32 23.532m0 0-43.701-46.833 76.772 1.794-33.07 45.04m6.972 37.496-25.947-13.782 18.974-23.715 6.973 37.497m-25.947-13.782-.728-22.72 19.702-.995-18.975 23.715M279.82 81.002l-16.265-4.512L326.93 93.48l-47.11-12.477m23.457 13.27-23.457-13.27 47.11 12.477-23.653.792m-34.304-2.676 10.847-10.593 23.457 13.27-34.304-2.677m-4.946-3.287-.472-11.818 16.265 4.512-15.793 7.306m4.946 3.287-4.946-3.287 15.793-7.306-10.847 10.593m43.132 171.247-24-47.828 43.702 46.833-19.702.995m-47.696-52.392-28.231-65.938 61.074 41.284-32.843 24.654m0 0-23.938.357-4.293-66.295 28.231 65.938m23.697 4.564-23.697-4.564 32.843-24.654-9.146 29.218m-38.814 51.773 4.797-52.572 10.32-3.765-15.117 56.337m0 0 15.117-56.337 23.697 4.564-38.814 51.773m4.797-52.572-13.618-3.408 23.938-.357-10.32 3.765m-23.925 68.442 19.128-15.87 17.863 8.439-36.99 7.43m32.882 36.995-32.883-36.994 36.99-7.431-4.107 44.425m-32.883-36.994-15.698-2.611 34.826-13.26-19.128 15.87m-31.377 35.553 31.377-35.552 32.883 36.994-64.26-1.442m0 0 15.679-38.163 15.698 2.61-31.377 35.553m193.418-111.67-10.211-23.527-7.135-41.838 17.346 65.365m-10.211-23.527-84.742 2.784 77.607-44.622 7.135 41.838m-6.77 86.896 7.903-31.593 9.078-31.776-16.98 63.37m-113.068 65.28 8.985-11.077 43.23-9.463-52.215 20.54m8.985-11.077-8.095-4.46 42.292-18.61-34.197 23.07m0 0 34.197-23.07 9.033 13.607-43.23 9.463m-8.985 11.076.89-15.537 8.095 4.46-8.985 11.077m65.364-17.346 11.26-18.498 36.444-29.436-47.704 47.934m11.26-18.498-6.974-37.497 43.418 8.061-36.445 29.436M236.178 144.512l-35.13-4.347 49.536-28.9-14.406 33.247m0 0 14.406-33.247 46.668 74.53-61.074-41.283m69.16 156.53-38.183-25.816 45.677 10.336-7.493 15.48m-38.184-25.816-17.863-8.44 38.814-51.772-20.951 60.212m44.95-12.384-44.95 12.384 20.95-60.212 24 47.828m.727 22.72-45.677-10.336 44.95-12.384.727 22.72m-15.58-99.766 50.865-66.18 26.742 21.558-77.607 44.622M150.515 253.46l-17.003-48.311 49.669 21.75-32.666 26.56m37.52-12.274-37.52 12.275 32.666-26.562 4.854 14.287m-37.182 29.32-.338-17.045 37.52-12.275-37.182 29.32m-.338-17.045.338 17.044-17.341-65.355 17.003 48.31m41.534 46.92-41.196-29.875 37.182-29.319 4.014 59.194m70.108 34.81 52.215-20.54 13.149 3.193-65.364 17.346m.89-15.537 4.108-44.425 38.184 25.817-42.292 18.608m-12.463-208.386 18.389-19.67 34.304 2.676-52.693 16.994m-109.058 73.264-8.014 20.62 16.975-63.379-8.961 42.76m0 0 12.563-38.72 39.693 24.01-52.256 14.71m44.332 38.005-44.332-38.006 52.256-14.708-7.924 52.714m-52.346-17.386 8.014-20.62 44.332 38.006-52.346-17.386m8.014-20.62 8.96-42.759 3.603 4.04-12.563 38.72m206.591-64.913L326.93 93.479l47.929 47.695-26.742-21.557m-44.84-25.346 23.653-.792 21.187 26.138-44.84-25.346m-6.025 91.525 6.025-91.525 44.84 25.346-50.865 66.179m0 0-46.668-74.531 52.693-16.994-6.025 91.525m77.973 84.112-43.418-8.06 51.32-23.533-7.902 31.593m-111.198-181.6-22.18-1.456 21.708-10.362.472 11.818m-13.443 22.957 13.443-22.957 4.946 3.287-18.389 19.67m0 0-8.736-24.413 22.179 1.456-13.443 22.957m-10.113 99.542-46.689-40.986 42.396-25.309 4.293 66.295m-46.689-40.986 7.266-29.656 35.13 4.347-42.396 25.309m-7.924 52.714 7.924-52.714 46.689 40.986-54.613 11.728m63.434 44.252-8.821-55.98 13.618 3.408-4.797 52.572m0 0-61.257-25.602 52.436-30.378 8.821 55.98m-55.51-96.966-27.85-30.64 35.116.984-7.266 29.656m-5.747 71.364-2.177-18.65 54.613-11.728-52.436 30.378m26.431 38.861-26.431-38.86 61.257 25.6-34.826 13.26M154.089 145.81l11.842-6.629 27.851 30.64-39.693-24.01m0 0-3.602-4.04 15.444-2.59-11.842 6.63m108.068 189.377-63.37-16.98 64.26 1.443-.89 15.537m-78.976-108.29-49.67-21.75 52.347 17.387-2.677 4.363m4.854 14.287-4.854-14.287 2.677-4.363 2.177 18.65m4.014 59.194-4.014-59.194 26.431 38.861-22.417 20.333m6.738 17.83-6.738-17.83 22.417-20.333-15.68 38.163m0 0-47.933-47.705 41.196 29.875 6.738 17.83"),
                    ui.tags.span("Database", class_="tooltip")
                 ),
        ),
            ui.tags.li(
                {"class": "icon-content"},
                ui.tags.a(
                    {"href": "https://www.linkedin.com/in/anasfplima/", 'aria-label':'Github', "data-social": "linkedin"},
                    ui.tags.div(class_ = "filled"),
                    ui.tags.svg(
                        {'xmlns':"http://www.w3.org/2000/svg",
                        'width':"16",
                        'height':"16",
                        'fill':"currentColor",
                        'viewBox':"0 0 640 640"},
                        ui.tags.Tag('path', d="M512 96L127.9 96C110.3 96 96 110.5 96 128.3L96 511.7C96 529.5 110.3 544 127.9 544L512 544C529.6 544 544 529.5 544 511.7L544 128.3C544 110.5 529.6 96 512 96zM231.4 480L165 480L165 266.2L231.5 266.2L231.5 480L231.4 480zM198.2 160C219.5 160 236.7 177.2 236.7 198.5C236.7 219.8 219.5 237 198.2 237C176.9 237 159.7 219.8 159.7 198.5C159.7 177.2 176.9 160 198.2 160zM480.3 480L413.9 480L413.9 376C413.9 351.2 413.4 319.3 379.4 319.3C344.8 319.3 339.5 346.3 339.5 374.2L339.5 480L273.1 480L273.1 266.2L336.8 266.2L336.8 295.4L337.7 295.4C346.6 278.6 368.3 260.9 400.6 260.9C467.8 260.9 480.3 305.2 480.3 362.8L480.3 480z"),
                    ),
                    ui.tags.span("LinkedIn", class_="tooltip")
                 ),
            ),     
        ),
        )
    ),
    window_title="ViruScope",
    id='main'
)




def server(input, output, session):
    # Create a temp directory just for THIS session
    session_temp_dir = make_temp_dir(project_dir)


    # =============================================================================
    # UI - Logo Dark mode toggle
    # =============================================================================
    @output
    @render.ui    
    def toggle_logo():
        """ toggle btw black/white logo based on dark/light mode """
        current_mode = input.mode()
        logo_img_by_mode = 'logo light.svg' if current_mode=='dark' else 'logo full.svg'
        logo_tag = ui.tags.img(
            src    = logo_img_by_mode,
            width = '200px',
            style  = 'margin-right: 10px;'
        )
        return logo_tag
    
    # =============================================================================
    # Initializing / Reactive value storage
    # =============================================================================
    vc = viruscope()
    al = arolit()
    primers_data = reactive.Value(None)
    status_text = reactive.Value("")
    calc_complete = reactive.Value(False)
    save_to = session_temp_dir
    loci_data = reactive.Value(None)
    combinations = reactive.Value(None)
    doi_data = reactive.Value(None)
    scraped_primers_list = reactive.Value(None)
    validated_primers_list = reactive.Value(None)
    scores_data = reactive.Value(None)

    # =============================================================================
    # Homepage
    # =============================================================================
    @reactive.effect
    @reactive.event(input.arolitcard)
    async def tab_switch():
        ui.update_navs(id="main", selected='arolit') 

    @reactive.effect
    @reactive.event(input.isopcard)
    async def tab_switch():
        ui.update_navs(id='main',selected='isop')           

    @reactive.effect
    @reactive.event(input.comboscard)
    async def tab_switch():
        ui.update_navs(id='main',selected='primer_combinations')

    @reactive.effect
    @reactive.event(input.utilitiescard)
    async def tab_switch():
        ui.update_navs(id='main',selected='utilities')   

    # =============================================================================
    # AROLit -  Fetch Metadata Section
    # =============================================================================
    @output
    @render.text
    def nbib_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.fetch_articles)
    def nbib_fetch():
        calc_complete.set(False)

        try:

            if not input.search_terms():
                ui.notification_show(
                    "Please enter a search query.",
                    duration=3,
                    type="error"
                )

                return
            
            with ui.Progress(min=0, max=1) as p:
                # status_text.set("Fetching articles...")
                p.set(message="Initializing", detail="Preparing to fetch articles...")

                query = input.search_terms()
                output_file = os.path.join(session_temp_dir, "articles.nbib")

                p.set(0.5, message="Fetching", detail="Querying PubMed database...")

                result = al.fetch_pubmed_medline(query, output_file)

                p.set(1, message="Completed", detail="Articles successfully fetched.")
                status_text.set(result)
                calc_complete.set(True)

        except Exception as e:
            status_text.set(f"An error occurred: {e}")

    @output
    @render.download(filename='articles.nbib')
    def download_nbib():
        file_path = os.path.join(session_temp_dir, 'articles.nbib')
        return file_path
    
    @output
    @render.ui
    def article_count():
        count = 0
        with open(os.path.join(session_temp_dir, 'articles.nbib'), 'r', encoding='utf-8') as file:
            for line in file:
                if line.strip() == '':
                    count += 1
        if count is not None or count != 0:
            return ui.tags.span([
                ui.tags.b("Total articles found: "),
                f"{count}"
            ],
            style="display: inline-block; margin-bottom: 10px;"
            )
        return ""

    @output
    @render.ui
    @reactive.event(input.fetch_articles)
    def metadata():
        while not calc_complete.get():
            return None
        file_path = os.path.join(session_temp_dir, 'articles.nbib')
        if os.path.exists(file_path):
            return (
                article_count,
                ui.tags.br(),
                ui.download_button("download_nbib", "Download .nbib file"),
            )

        return "No articles found."
    # =============================================================================
    # AROLit -  Primer Scraping Section
    # =============================================================================
    @output
    @render.text
    def scrape_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.scrape_primers)
    def run_scrape():
        calc_complete.set(False)

        try:
            if not input.for_scraping():
                ui.notification_show(
                    "Please upload at least one PDF file for primer retrieval.",
                    duration=3,
                    type="error"
                )
                return

            with ui.Progress(min=0, max=1) as p:
                #status_text.set("Scraping primers...")
                p.set(0, message="Initializing", detail="Preparing to scrape primers...")

                pdf_files = [file["datapath"] for file in input.for_scraping()]
                names = [file["name"] for file in input.for_scraping()]

                p.set(0.5, message="Scraping", detail="Extracting primers from PDFs...")

                result = al.oligos_to_csv_mult(pdf_files, names, os.path.join(session_temp_dir, "primers.csv"))

                scraped_primers_list.set({
                    "primers": al.primers_dict,
                    "count": len(al.primers_dict)
                })

                if input.sql_db():
                    try:
                        status_text.set("Saving primers to SQL database...")
                        al.csv_to_sql(os.path.join(session_temp_dir, "primers.csv"), os.path.join(session_temp_dir, "primers.db"))
                    except Exception as e:
                        status_text.set(f"An error occurred while saving to SQL database: {e}")
                        pass

                p.set(1, message="Completed", detail="Primers successfully scraped.")
                status_text.set(result)
                calc_complete.set(True)

        except Exception as e:
            status_text.set(f"An error occurred: {e}")

    @output
    @render.download(filename='primers.csv')
    def download_arolit_primers():
        file_path = os.path.join(session_temp_dir, 'primers.csv')
        return file_path
    
    @output
    @render.download(filename='primers.db')
    def download_db():
        file_path = os.path.join(session_temp_dir, 'primers.db')
        return file_path

    @output
    @render.data_frame
    def arolit_table():
        df = pd.DataFrame(scraped_primers_list()["primers"], columns=["Sequence", "Source"])
        if df is None:
            return None
        else:
            df = df.head(1000)
            return render.DataTable(df,styles={"class": "datagrid"}, summary= 'Showing the first 1000 rows. Download the CSV for the complete primer list.')
        
    @output
    @render.ui
    @reactive.event(input.scrape_primers)
    def scraped_primers():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(scraped_primers_list()["primers"], columns=["Sequence", "Source"])
        if df is None:
            return None
        return (
            ui.tags.span([
                ui.tags.b("Total of primers found: "),
                f"{scraped_primers_list()['count']}"
            ],
            style="display: inline-block; margin-bottom: 10px;"
            ),
            arolit_table,
            ui.div(
                [
                ui.download_button("download_arolit_primers", "Download CSV"),
                ui.download_button("download_db", "Download SQL Database") if input.sql_db() else None
            ],
            style="display: flex; justify-content: space-between;",)
    )
        #     ui.download_button("download_arolit_primers", "Download Primers (CSV)"),
        #     ui.download_button("download_db", "Download SQL Database") if input.sql_db() else None
        # )

    # =============================================================================
    # AROLit - Primer Validation Section
    # =============================================================================

    @output
    @render.text
    def validate_status():
        return status_text.get()

    @reactive.effect
    @reactive.event(input.validate_primers)
    def run_validation():
        calc_complete.set(False)

        # Input checks
        if not input.for_validation():
            ui.notification_show(
                "Please upload a primer database for validation.",
                duration=3,
                type="error"
            )
            return
        if not (input.validate_ref() or input.validate_align()):
            ui.notification_show(
                "BLAST selected but no genome file provided.",
                duration=3,
                type="error"
            )
            return
        
        try:
            with ui.Progress(min=0, max=1) as p:
                # status_text.set("Validating primers...")
                p.set(0, message="Initializing", detail="Preparing to validate primers...")

                # Prepare inputs
                primer_file = input.for_validation()[0]["datapath"]
                print( f"Validating primers from: {primer_file}" )
                method = input.blast_method()
                print( f"Using validation method: {method}" )
                if input.blast_method() == 'genome_ref':
                    genome_path = input.validate_ref()[0]["datapath"]
                elif input.blast_method() == 'genome_aligned':
                    genome_path = input.validate_align()[0]["datapath"]
                print( f"Using genome file: {genome_path}" )
                save_to = session_temp_dir
                
                p.set(0.3, message="Validating", detail="Making BLAST database...")

                if method == 'genome_ref':
                    result1 = vc.make_blastdb_ref(save_to, genome_path)
                    db_path = os.path.join(save_to, "blastref")
                elif method == 'genome_aligned':
                    result1 = vc.make_blastdb_alig(save_to, genome_path)
                    db_path = os.path.join(save_to, "blastalign")

                status_text.set(result1)
                p.set(0.6, message="Validating", detail="Running BLAST...")

                vc.export_for_score(
                    input_db=primer_file,
                    output_csv=os.path.join(save_to, "primers_for_scoring.csv"),
                )

                vc.primers_to_fasta(
                    primers_csv=os.path.join(save_to, "primers_for_scoring.csv"),
                    fasta_output=os.path.join(save_to, "primers_for_scoring.fasta"),
                )

                result2 = vc.run_blast(
                    primers_fasta=os.path.join(save_to, "primers_for_scoring.fasta"),
                    db_path=db_path,
                    save_to=save_to,
                    )
                
                status_text.set(result2)
                p.set(0.8, message="Parsing", detail="Parsing BLAST results...")
                

                result3 = vc.parse_blast(
                    results_file=os.path.join(save_to, "primer_hits.txt"),
                    input_csv=os.path.join(save_to, "primers_for_scoring.csv"),
                    output_csv=os.path.join(save_to, "validated_primers.csv"),
                )
                
                vc.GC_MT_calculations('validated_primers')
                vc.primer_type('validated_primers')

                validated_primers_list.set({
                    "primers": vc.validated_primers,
                    "count": len(vc.validated_primers)
                })


                p.set(1.0, message="Completed", detail="Primer validation completed successfully.")
                calc_complete.set(True)

        except Exception as e:
            status_text.set(f"An error occurred: {e}")

    # Output table
    @output
    @render.data_frame
    def validated_table():
        df = pd.DataFrame(validated_primers_list()["primers"])
        if df is None:
            return None
        else:
            df = df.head(1000)
            return render.DataTable(df,styles={"class": "datagrid"}, summary= 'Showing the first 1000 rows. Download the CSV for the complete primer list.')
        
    @output
    @render.download(filename='validated_primers.csv')
    def download_validated_primers():
        df = pd.DataFrame(validated_primers_list()["primers"])
        df.to_csv(os.path.join(session_temp_dir, 'validated_primers.csv'), index=False)
        return os.path.join(session_temp_dir, 'validated_primers.csv')
    
    @output
    @render.download(filename='validated_primers.fasta')
    def download_validated_fasta():
        utilities.export_primers_to_fasta(validated_primers_list()["primers"], os.path.join(session_temp_dir, 'validated_primers.fasta'))
        file_path = os.path.join(session_temp_dir, 'validated_primers.fasta')
        return file_path

    @output
    @render.ui
    @reactive.event(input.validate_primers)
    def valid_primers():
        while not calc_complete.get():
            return None
        df = pd.read_csv(os.path.join(session_temp_dir, "validated_primers.csv"))
        if df is None:
            return None
        return (
            ui.tags.span([
                ui.tags.b("Total of primers found: "),
                f"{validated_primers_list()['count']}"
            ],
            style="display: inline-block; margin-bottom: 10px;"
            ),
            validated_table,
            ui.div(
                [
                ui.download_button("download_validated_primers", "Download CSV"),
                ui.download_button("download_validated_fasta", "Download FASTA")
            ],
            style="display: flex; justify-content: space-between;",)
    )
        
        
    # =============================================================================
    # ISOP - Generate Primers Section
    # =============================================================================

    @output
    @render.text
    def blast_status():
        return status_text.get()

    @reactive.effect
    @reactive.event(input.generate_primers)
    def run_generate():
        calc_complete.set(False)

        # Input checks
        if not input.genome_file():
            ui.notification_show(
                "Please upload a genome file.",
                duration=3,
                type="error"
            )
            return
        if input.blast_filter() and not input.align_ref():
            ui.notification_show(
                "BLAST selected but no reference genome provided.",
                duration=3,
                type="error"
            )
            return

        try:
            with ui.Progress(min=0, max=1) as p:
                # status_text.set("Generating primers...")
                p.set(0, message="Initializing", detail="Preparing to generate primers...")

                # Prepare inputs
                genome_path = input.genome_file()[0]["datapath"]
                name = input.primer_name()
                r_start = input.min_length()
                r_end = input.max_length()
                step = 1
                blast = input.blast_filter()
                save_to = session_temp_dir


                if blast:
                    alignment_path = input.align_ref()[0]["datapath"]

                    p.set(0.3, message="Generating", detail="Running primer generation...")

                    db_status = vc.make_blastdb_alig(save_to, alignment_path)

                    p.set(0.4, message="Generating", detail="Running BLAST...")
                    
                    result = vc.generate_primers(
                        genome_file=genome_path,
                        name=name,
                        save_to=save_to,
                        output_file="primer_output.txt",
                        range_start=r_start,
                        range_end=r_end,
                        step=step,
                        blast=True,
                        db_path=os.path.join(save_to, "blastalign")
                    )
                    status_text.set(f"{db_status}\n\n{result}")

                else:
                    p.set(0.3, message="Generating", detail="Running primer generation...")
                    result = vc.generate_primers(
                        genome_file=genome_path,
                        name=name,
                        save_to=save_to,
                        output_file="primer_output.txt",
                        range_start=r_start,
                        range_end=r_end,
                        step=step,
                        blast=False
                    )
                    status_text.set(result)


                try:
                    p.set(0.6, message="Calculating", detail="Running GC/MT calculations...")

                    status_text.set("Calculating GC content and melting temperatures...")

                    vc.GC_MT_calculations('silico_primers' if not blast else 'filtered_silico_primers')

                    status_text.set(f"GC and melting temperature calculations completed.")

                except Exception as e:
                    status_text.set(f"An error occurred during GC/MT calculations: {e}")

                # Store primers in reactive
                primers_data.set({
                    "before": vc.silico_primers,
                    "after": vc.filtered_silico_primers,
                    "blast_ran": blast,
                    "count": len(vc.filtered_silico_primers if blast else vc.silico_primers)
                })

                p.set(1.0, message="Completed", detail="Primer generation completed successfully.")
                calc_complete.set(True)
                status_text.set("Primer generation completed successfully.")
            
        except Exception as e:
            status_text.set(f"An error occurred: {e}")

    # Output table
    @output
    @render.data_frame
    def primer_outputs_table():
        if input.blast_filter():
            df = pd.DataFrame(primers_data()["after"])
        else:
            df = pd.DataFrame(primers_data()["before"])
        if df is None:
            return None
        else:
            df = df.head(1000)
            df = df.rename(columns={
                'GC Content': 'GC Content (%)',
                'Melting Temperature': 'Melting Temperature (°C)'
                })
            # df = df.drop(columns=['PPI', 'PPI3', 'Score'])
        return render.DataTable(df,styles={"class": "datagrid"}, summary= 'Showing the first 1000 rows. Download the CSV for the complete primer list.')

    # Output count
    @output
    @render.ui
    def primer_count():
        count = primers_data()["count"]
        if count is not None:
            return ui.tags.span([
                ui.tags.b("Total primers: "),
                f"{count}"
            ],
            style="display: inline-block; margin-bottom: 10px;"
            )
        return ""
    
    # CSV download
    @output
    @render.download
    def download_csv():
        if input.blast_filter():
            filename = "primersafterblast.csv"
            df = pd.DataFrame(primers_data()["after"])
            # df = df.drop(columns=['PPI', 'PPI3', 'Score'])
            df.to_csv(os.path.join(session_temp_dir, filename), index=False)
            return os.path.join(session_temp_dir, filename)
        else:
            filename = "primersbeforeblast.csv"
            df = pd.DataFrame(primers_data()["before"])
            # df = df.drop(columns=['PPI', 'PPI3', 'Score'])
            df.to_csv(os.path.join(session_temp_dir, filename), index=False)
            return os.path.join(session_temp_dir, filename)

    # FASTA download
    @output
    @render.download
    def download_fasta():
        if input.blast_filter():
            fasta_file = os.path.join(session_temp_dir, "primersafterblast.fasta")
        else:
            fasta_file = os.path.join(session_temp_dir, "primersbeforeblast.fasta")
        return fasta_file

    # Dynamic UI area
    @output
    @render.ui
    @reactive.event(input.generate_primers)
    def primer_outputs():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(primers_data()["before"] if not input.blast_filter() else primers_data()["after"])
        if df is None:
            return ui.div()
        return ui.div(
            primer_count,
            primer_outputs_table,
            ui.download_button("download_csv", "Download Primers (CSV)"),
            ui.download_button("download_fasta", "Download Primers (FASTA)")
        )
    
    # =============================================================================
    # ISOP - Conservation Scores Section
    # =============================================================================
    @output
    @render.text
    def CS_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.calculate_CS)
    def run_CS():
        calc_complete.set(False)

        try:

            if not input.align_file():
                ui.notification_show(
                    "Please upload an alignment file.",
                    duration=3,
                    type="error"
                )
                return
            if input.upload_method() == "upload" and not input.primer_csv():
                ui.notification_show(
                    "Please upload a primer CSV file.",
                    duration=3,
                    type="error"
                )
                return
            
            if input.upload_method() == "current" and primers_data() is None:
                ui.notification_show(
                    "No primers available for conservation score calculations.",
                    duration=3,
                    type="error"
                )
                return

            with ui.Progress(min=0, max=1) as p:

                # Input checks
                align_path = input.align_file()[0]["datapath"]

                if input.upload_method() == "upload":
                    
                    p.set(0, message="Loading", detail="Loading uploaded primers...")

                    primers_path = input.primer_csv()[0]["datapath"]
                    vc.load_primers_from_csv(primers_path)
                    primers_df = 'uploaded_primers'
                else:
                    p.set(0, message="Initializing", detail="Preparing to calculate conservation scores...")
                    if primers_data()['blast_ran'] == False:
                        primers_df = 'silico_primers'
                    else:
                        primers_df = 'filtered_silico_primers'
                
                if primers_df is None:
                    ui.notification_show(
                        "No ",ui.tags.i("in silico"), " primers available for conservation score calculations.",
                        duration=3,
                        type="error"
                    )
                    return

                p.set(0.4, message="Calculating", detail="This may take a while...")

                # Calculate conservation scores
                cs_results = vc.parse_primers_app(align_path, primers_df, session_temp_dir, status_callback=status_text.set)

                # Store results in reactive without overwriting existing data
                current = primers_data() or {}
                current["conservation_scores"] = cs_results
                scores_data.set(current)

                p.set(1, message="Complete", detail="Conservation scores calculated successfully.")
                calc_complete.set(True)
                status_text.set("Conservation scores calculated successfully.")

        except Exception as e:
            status_text.set(f"An error occurred: {e}")

    # Output table
    @output
    @render.data_frame
    def conservation_scores_table():
        file_path = os.path.join(session_temp_dir, 'conservation_scores.csv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            if df.empty:
                return None
            else:
                df = df.head(1000)
                df = df[['Primer ID', 'Sequence', 'PPI', 'PPI3', 'Score']]
                return render.DataTable(df, styles={"class": "datagrid"}, summary='Showing the first 1000 rows. Download the CSV for the complete primer list.')

    # CSV download
    @output
    @render.download(filename='conservation_scores.csv')
    def download_CS():
        file_path = os.path.join(session_temp_dir, 'conservation_scores.csv')
        return file_path

    # Output UI
    @output
    @render.ui
    @reactive.event(input.calculate_CS)
    def CS_output():
        while not calc_complete.get():
            return None
        df = pd.read_csv(os.path.join(session_temp_dir, 'conservation_scores.csv'))
        if df.empty:
            return ui.div("No conservation scores available.")
        else:
            return ui.div(
                conservation_scores_table,
                ui.download_button("download_CS", "Download Conservation Scores (CSV)")
            )
        
    # =============================================================================
    # Gene extraction and loci mapping Section
    # =============================================================================
    @output
    @render.text
    def loci_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.retrieve_loci)
    def loci_retrieval():
        calc_complete.set(False)
        loci_info_mapped = None

        try:
            # Input checks
            if input.retrieve_method() == "auto" and not input.email():
                ui.notification_show(
                    "Please enter your email address.",
                    duration=3,
                    type="error"
                )
                return
            if input.retrieve_method() == "auto" and not input.ncbi_id():
                ui.notification_show(
                    "Please enter an NCBI ID.",
                    duration=3,
                    type="error"
                )
                return
            if input.retrieve_method() == "upload" and not input.genbank_file():
                ui.notification_show(
                    "Please upload a GenBank file.",
                    duration=3,
                    type="error"
                )
                return
            if input.map_to() and not input.alignment_file():
                ui.notification_show(
                    "Please upload an alignment file.",
                    duration=3,
                    type="error"
                )
                return
            if input.map_to() and not input.reference_file():
                ui.notification_show(
                    'Please upload a FASTA file with the reference genome.',
                    duration=3,
                    type="error"
                )
                return

            status_text.set("Retrieving gene loci...")

            if input.retrieve_method() == "upload":             
                genbank_path = input.genbank_file()[0]["datapath"]
                ncbi_id = None
                email = None

            else:
                genbank_path = None
                ncbi_id = input.ncbi_id()
                email = input.email()

            # Retrieve loci
            loci_info = vc.extract_loci_from_genbank(
                genbank_file=genbank_path,
                ncbi_id=ncbi_id,
                email=email
            )
            print(f"Extracted loci: {loci_info}")  # Debug

            # If user wants to map to alignment
            if input.map_to():
                status_text.set("Mapping loci to alignment...")

                alignment_file = input.alignment_file()[0]["datapath"]
                reference_file = input.reference_file()[0]["datapath"]

                loci_info_mapped = vc.map_ref_coords_to_alignment(
                    reference_file,  # unaligned reference
                    alignment_file,  # aligned reference
                    loci_info
                )
                print(f"Mapped loci: {loci_info_mapped}")  # Debug

            loci_data.set({
                'ref': loci_info,
                'mapped': loci_info_mapped
            })

            calc_complete.set(True)
            status_text.set("Gene loci retrieval completed successfully.")

        except Exception as e:
            import traceback
            traceback.print_exc()
            status_text.set(f"An error occurred: {e}")
    
    # Output table (reference)
    @output
    @render.data_frame
    def ref_loci_table():
        df = pd.DataFrame(loci_data()['ref'])
        df = df.rename(columns={'name':'Gene', 'aligned_start':'Start', 'aligned_end':'End'})
        if df.empty:
            return ui.div("No reference loci available.")
        return render.DataTable(df, styles={"class": "datagrid"}, height="auto")

    # Output table (mapped)
    @output
    @render.data_frame
    def mapped_loci_table():
        df = pd.DataFrame(loci_data()['mapped'])
        df = df.rename(columns={'name':'Gene', 'start':'Start', 'end':'End'})
        if df.empty:
            return ui.div("No mapped loci available.")
        return render.DataTable(df, styles={"class": "datagrid"}, height="auto")
    
    # Download buttons
    @output
    @render.download()
    def download_ref_loci():
        df = pd.DataFrame(loci_data()['ref'])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'loci_reference.csv'), index=False)
        return os.path.join(session_temp_dir, 'loci_reference.csv')

    @output
    @render.download()
    def download_mapped_loci():
        df = pd.DataFrame(loci_data()['mapped'])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'loci_mapped.csv'), index=False)
        return os.path.join(session_temp_dir, 'loci_mapped.csv')
    
    # Output UI
    @output
    @render.ui
    @reactive.event(input.retrieve_loci)
    def lociref_output():
        while not calc_complete.get():
            return None
        return ui.card(
            ui.card_header("Loci - Reference genome"),
            ui.card_body(ref_loci_table),
            ui.download_button("download_ref_loci", "Download Reference Loci (CSV)")
        )

    @output
    @render.ui
    @reactive.event(input.retrieve_loci)
    def locialign_output():
        if input.map_to():
            while not calc_complete.get():
                return None
            return ui.card(
                ui.card_header("Loci - Mapped genome"),
                ui.card_body(mapped_loci_table),
                ui.download_button("download_mapped_loci", "Download Mapped Loci (CSV)")
            )
        
    # =============================================================================
    # Annotate primers Section
    # =============================================================================
    @reactive.effect
    @reactive.event(input.cs_switch)
    async def tab_switch():
        # ui.update_navs(id='main',selected='isop')
        # await asyncio.sleep(0.05)
        ui.update_navs(id='combinations_tabs',selected='conservation_scores')

    @output
    @render.text
    def annotated_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.annotate_primers)
    def annotate():
        calc_complete.set(False)
        try:
            # Input checks
            if not input.loci_file():
                ui.notification_show(
                    "Please upload a loci file.",
                    duration=3,
                    type="error"
                )
                return
            if not input.primer_file():
                ui.notification_show(
                    "Please upload a primer file.",
                    duration=3,
                    type="error"
                )
                return

            status_text.set("Annotating primers with loci information...")

            loci_path = input.loci_file()[0]["datapath"]
            primer_path = input.primer_file()[0]["datapath"]
            loci_list = vc.load_primers_from_csv(loci_path)
            vc.load_primers_from_csv(primer_path)
            primer_list = 'uploaded_primers'

            # Annotate primers
            annotated_primers = vc.annotate_primers_with_locus('uploaded_primers', loci_list)

            # Store results in reactive
            current = primers_data() or {}
            current["annotated_primers"] = annotated_primers
            primers_data.set(current)

            calc_complete.set(True)
            status_text.set("Primer annotation completed successfully.")

        except Exception as e:
            status_text.set(f"An error occurred: {e}")
    
    @output
    @render.data_frame
    def annotated_primers_table():
        df = pd.DataFrame(primers_data()["annotated_primers"])
        if df.empty:
            return ui.div("No primers to annotate.")
        else:
            df = df.head(1000)
            df = df[['Primer ID', 'Sequence', 'Locus']]
            return render.DataTable(df, styles={"class": "datagrid"}, summary='Showing the first 1000 rows. Download the CSV for the complete primer list.')

    @output
    @render.download()
    def download_annotated_primers():
        df = pd.DataFrame(primers_data()["annotated_primers"])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'primers_annotated.csv'), index=False)
        return os.path.join(session_temp_dir, 'primers_annotated.csv')
    
    @output
    @render.ui
    @reactive.event(input.annotate_primers)
    def annotated_primers():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(primers_data()["annotated_primers"])
        if df.empty:
            return ui.div("No annotated primers available.")
        return ui.div(
            annotated_primers_table,
            ui.download_button("download_annotated_primers", "Download Annotated Primers (CSV)")
        )
    
    # =============================================================================
    # Primer Combination section
    # =============================================================================
    @reactive.effect
    @reactive.event(input.annotation_switch)
    async def tab_switch():
        ui.update_navs(id='main',selected='primer_combinations')
        await asyncio.sleep(0.05)
        ui.update_navs(id='combinations_tabs',selected='annotate_primers')

    
    @output
    @render.text
    def combinations_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.explore_combinations)
    async def calc_combs():
        calc_complete.set(False)

        try:
            # Input checks
            if not input.primers_ready():
                ui.notification_show(
                    "Please upload a CSV file containing primers.",
                    duration=3,
                    type="error"
                )
                return

            with ui.Progress(min=0, max=1) as p:
                # status_text.set("Calculating primer combinations...")
                p.set(0, message="Initializing", detail="Preparing to calculate primer combinations...")

                primers_path = input.primers_ready()[0]["datapath"]
                custom_params = input.custom_parameters()
                amp_min = input.amp_min() if custom_params else 100
                amp_max = input.amp_max() if custom_params else 1000
                gc_min = input.gc_min() if custom_params else 40
                gc_max = input.gc_max() if custom_params else 60
                mt_min = input.mt_min() if custom_params else 55
                mt_max = input.mt_max() if custom_params else 65
                min_cs = input.min_cs() if custom_params else 99

                p.set(0.3, message="Loading", detail="Loading uploaded primers...")

                # Load primers
                uploaded_primers = vc.load_primers_from_csv(primers_path)
                

                p.set(0.5, message="Processing", detail="Calculating primer combinations...")

                # Calculate combinations
                vc.calculate_primer_combinations(
                    'uploaded_primers',
                    min_amplicon_length=amp_min,
                    max_amplicon_length=amp_max,
                    conservation_score_threshold=min_cs,
                    min_GC=gc_min,
                    max_GC=gc_max,
                    min_tm=mt_min,
                    max_tm=mt_max,
                )

                # Store results in reactive
                combinations.set({
                    'combos': vc.primer_combinations,
                    'count': len(vc.primer_combinations)
                })


                p.set(0.8, message="Calculating", detail="Obtaining folding scores for primer combinations...")

                # Fold scores
                vc.fold_scores_for_combinations('primer_combinations')
                current = combinations.get() or {}
                current["combos"] = vc.primer_combinations
                combinations.set(current)

                p.set(1, message="Complete", detail="Primer combinations calculated successfully.")
                calc_complete.set(True)
                status_text.set("Primer combinations calculated successfully.")

        except Exception as e:
            status_text.set(f"An error occurred: {e}")
    
    @output
    @render.data_frame
    def comb_table():
        df = pd.DataFrame(combinations()["combos"])
        if df.empty:
            return ui.div("No primer combinations available.")
        else:
            df = df.head(1000)
            df = df[['Primer ID_fwd', 'Primer ID_rev', 'Amplicon_length','GC_content', 'Melting Temperature', 'Combination_score','No-fold score']]
            return render.DataTable(df, styles={"class": "datagrid"}, summary='Showing the first 1000 rows. Download the CSV for the complete primer list.')

    @output
    @render.download()
    def download_combos():
        df = pd.DataFrame(combinations()["combos"])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'primers_combinations.csv'), index=False)
        return os.path.join(session_temp_dir, 'primers_combinations.csv')

    # Output count
    @output
    @render.ui
    def pair_count():
        count = combinations()["count"]
        if count is not None:
            return ui.tags.span([
                ui.tags.b("Total primer combinations: "),
                f"{count}"
            ],
            style="display: inline-block; margin-bottom: 10px;"
            )
        return ""
    
    @output
    @render.ui
    @reactive.event(input.explore_combinations)
    def primer_combinations():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(combinations()["combos"])
        if df.empty:
            return ui.div("No primer combinations available.")
        return ui.div(
            pair_count,
            comb_table,
            ui.download_button("download_combos", "Download Primer Combinations (CSV)")
        )
    
    @render.plot(alt="GC histogram")
    @reactive.event(input.explore_combinations)  
    def gc_histogram():
        while not calc_complete.get():
            return None
        if input.show_histograms():
            df = pd.DataFrame(combinations()["combos"])
            gc = df["GC_content"]

            fig, ax = plt.subplots()
            ax.hist(gc, bins=20, color='skyblue', edgecolor='black')
            ax.set_title("GC Content Distribution of all primers")
            ax.set_xlabel("GC Content (%)")
            ax.set_ylabel("Number of primers")
            ax.grid(True)

            return fig
        
        else:
            return None
    
    @render.plot(alt="TM histogram")
    @reactive.event(input.explore_combinations)
    def tm_histogram():  
        while not calc_complete.get():
            return None
        if input.show_histograms():
            df = pd.DataFrame(combinations()["combos"])
            tm = df["Melting Temperature"]

            fig, ax = plt.subplots()
            ax.hist(tm, bins=20, color='skyblue', edgecolor='black')
            ax.set_title("Melting Temperature Distribution of all primers")
            ax.set_xlabel("Melting Temperature (°C)")
            ax.set_ylabel("Number of primers")
            ax.grid(True)

            return fig
        
        else:
            return None
    
    
    # =============================================================================
    # Utilities - Fetch DOIs and PMIDs
    # =============================================================================
    @output
    @render.text
    def fetch_status():
        return status_text.get()
    
    @reactive.effect
    @reactive.event(input.fetch_dois)
    def fetch_():
        calc_complete.set(False)


        if not input.nbib_file():
            ui.notification_show(
                "Please upload a .nbib file.",
                duration=3,
                type="error"
            )
            return
        
        else:
            nbib_file = input.nbib_file()[0]["datapath"]
            status_text.set("Fetching DOIs and PMIDs...")
            try:
                doi_dict, pmid_dict = utilities.fetch_doi_dict(nbib_file)
                doi_data.set({
                    "DOIs": [key for key in doi_dict],
                    "PMIDs": [key for key in pmid_dict]
                })

                calc_complete.set(True)
                status_text.set('Retrieved DOIs and/or PMIDs successfully.')

            except Exception as e:
                status_text.set(f"An error occurred while fetching DOIs/PMIDs: {e}")

    @output
    @render.data_frame
    def doi_table():
        df = pd.DataFrame(doi_data()['DOIs'], columns=["DOIs"])
        if df.empty:
            return ui.div("No DOIs found in the .nbib file.")
        df = df.head(1000)
        return render.DataTable(df, styles={"class": "datagrid"}, summary='Showing the first 1000 rows. Download the CSV for the complete list.')

    @output
    @render.data_frame
    def pmid_table():
        df = pd.DataFrame(doi_data()['PMIDs'], columns=["PMIDs"])
        if df.empty:
            return ui.div("No PMIDs found in the .nbib file.")
        df = df.head(1000)
        return render.DataTable(df, styles={"class": "datagrid"}, summary='Showing the first 1000 rows. Download the CSV for the complete list.')

    @output
    @render.download()
    def download_dois():
        df = pd.DataFrame(doi_data()['DOIs'], columns=["DOIs"])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'dois_list.csv'), index=False)
        return os.path.join(session_temp_dir, 'dois_list.csv')
    
    @output
    @render.download()
    def download_pmids():
        df = pd.DataFrame(doi_data()['PMIDs'], columns=["PMIDs"])
        if df.empty:
            return None
        df.to_csv(os.path.join(session_temp_dir, 'pmids_list.csv'), index=False)
        return os.path.join(session_temp_dir, 'pmids_list.csv')

    @output
    @render.ui
    @reactive.event(input.fetch_dois)
    def doi_output():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(doi_data()['DOIs'], columns=["DOIs"])
        if df.empty:
            return ''
        else:
            return ui.div(
                doi_table,
                ui.download_button("download_dois", "Download DOIs (CSV)")
            )

    @output
    @render.ui
    @reactive.event(input.fetch_dois)
    def pmid_output():
        while not calc_complete.get():
            return None
        df = pd.DataFrame(doi_data()['PMIDs'], columns=["PMIDs"])
        if df.empty:
            return ''
        else:
            return ui.div(
                pmid_table,
                ui.download_button("download_pmids", "Download PMIDs (CSV)")
            )

    # =============================================================================
    # Utilities - Count sequences
    # =============================================================================

    @output
    @render.text
    @reactive.event(input.count_sequences)
    def sequence_count():
        if not input.fasta_file():
            ui.notification_show(
                "Please upload a FASTA file.",
                duration=3,
                type="error"
            )
            return
        else:
            fasta_file = input.fasta_file()[0]["datapath"]
            try:
                count = utilities.count_sequences(fasta_file)
                return count
            except Exception as e:
                return f"\nAn error occurred while counting sequences: {e}" 

    # =============================================================================
    # Utilities - Reorder sequences
    # =============================================================================
    @output
    @render.ui
    def download_reordered_ui():
        # Check if inputs are ready
        result = None
        is_ready = bool(input.reference_id()!="") and bool(input.fasta_file_reorder()!=None)
        if is_ready:
            fasta_file = input.fasta_file_reorder()[0]["datapath"]
            ref_id = input.reference_id()
            result = utilities.reorder_sequences(fasta_file,'reordered_sequences.fasta', ref_id, session_temp_dir)
        
        # Apply CSS to disable button when not ready
        style = "" if result!=None else "pointer-events: none; opacity: 0.5; cursor: not-allowed;"
        
        # Tooltip text if disabled
        tooltip = "No sequence found. Download unavailable." if is_ready and (result == None) else ""

        return ui.tags.div(
            ui.download_button(
                "download_reordered",
                "Download FASTA",
                width="fit-content",
                style=style
            ),
            title=tooltip
        )
    
    @output
    @render.download()
    def download_reordered():
        return os.path.join(session_temp_dir, 'reordered_sequences.fasta')


    # =============================================================================
    # Utilities - Remove duplicates
    # =============================================================================
    @output
    @render.ui
    def download_nodupes_ui():
        # Check if inputs are ready
        result = None
        is_ready = bool(input.fasta_file_dupes()!=None)
        if is_ready:
            fasta_file = input.fasta_file_dupes()[0]["datapath"]
            result = utilities.remove_dupes(fasta_file, 'noduplicates.fasta', session_temp_dir)
        
        # Apply CSS to disable button when not ready
        style = "" if result!=None else "pointer-events: none; opacity: 0.5; cursor: not-allowed;"
        
        # Tooltip text if disabled
        tooltip = "No file uploaded. Download unavailable." if not is_ready else ""

        return ui.tags.div(
            ui.download_button(
                "download_nodupes",
                "Download FASTA",
                width="fit-content",
                style=style
            ),
            title=tooltip
        )
    
    @output
    @render.download()
    def download_nodupes():
        return os.path.join(session_temp_dir, 'noduplicates.fasta')

    # =============================================================================
    # Utilities - Extract sequence
    # =============================================================================
    @output
    @render.ui
    def download_extract_ui():
        # Check if inputs are ready
        result = None
        is_ready = bool(input.sequence_id()!="") and bool(input.fasta_file_extract()!=None)
        if is_ready:
            fasta_file = input.fasta_file_extract()[0]["datapath"]
            sequence_id = input.sequence_id()
            output_file = f'{sequence_id}.fasta'
            result = utilities.extract_sequence(fasta_file, output_file, sequence_id, session_temp_dir)
        
        # Apply CSS to disable button when not ready
        style = "" if result!=None else "pointer-events: none; opacity: 0.5; cursor: not-allowed;"
        
        # Tooltip text if disabled
        tooltip = "No sequence found. Download unavailable." if is_ready and (result == None) else ""

        return ui.tags.div(
            ui.download_button(
                "download_extract",
                "Download FASTA",
                width="fit-content",
                style=style
            ),
            title=tooltip
        )
    
    @output
    @render.download()
    def download_extract():
        sequence_id = input.sequence_id()
        output_file = f'{sequence_id}.fasta'
        return os.path.join(session_temp_dir, output_file)
        # try:
        #     fasta_file = input.fasta_file_extract()[0]["datapath"]
        #     sequence_id = input.sequence_id()
        #     output_file = f'{sequence_id}.fasta'

        #     result = utilities.extract_sequence(fasta_file, output_file, sequence_id, session_temp_dir)
        #     if result is None:
        #         raise Exception(f"Sequence ID '{sequence_id}' not found in the provided FASTA file.")
        #     else:
        #         return os.path.join(session_temp_dir, output_file)

        # except Exception as e:
        #     raise Exception(f"An error occurred while extracting sequence: {e}")

    # =============================================================================
    # Utilities - Venn Diagram
    # =============================================================================

    @output
    @render.plot
    @reactive.event(input.generate_venn)
    def venn_diagram():
        if not input.fasta_file_venn1() or not input.fasta_file_venn2():
            ui.notification_show(
                "Please upload two FASTA files with primer sequences.",
                duration=3,
                type="error"
            )
            return None

        fasta_file1 = {str(record.seq) for record in SeqIO.parse(input.fasta_file_venn1()[0]["datapath"], "fasta")}
        fasta_file2 = {str(record.seq) for record in SeqIO.parse(input.fasta_file_venn2()[0]["datapath"], "fasta")}
        fasta_label1 = input.fasta_file_venn1()[0]["name"]
        fasta_label2 = input.fasta_file_venn2()[0]["name"]

        if fasta_file1 and fasta_file2:
            fig, ax = plt.subplots()
            venn2(
                [fasta_file1, fasta_file2],
                set_labels=(fasta_label1, fasta_label2),
                set_colors=('#007bc2', '#96ec85'),
                ax=ax)
            ax.set_title("Primer Overlap")
            return fig            
        else:
            return None
        
    # =============================================================================
    # SESSION END CLEANUP
    # =============================================================================
    session.on_ended(lambda: os.chdir(project_dir))
    session.on_ended(lambda: safe_delete_folder(session_temp_dir))


app = App(app_ui, server, static_assets=www_dir)
# atexit.register(shutil.rmtree, last_temp_dir)