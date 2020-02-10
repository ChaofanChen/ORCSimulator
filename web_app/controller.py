# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:01:36 2019

@author: Chaofan
"""
from model import InputForm
from flask import Flask, render_template, request
from compute import compute

app = Flask(__name__)

@app.route('/chaofan', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
#
    if request.method == 'POST' and form.validate():
        working_type = form.fluid_type.data
        result = compute(working_type, form.Q_w.data, form.Q_s.data,
                         form.T_b_p.data, form.T_b_i.data,
                         form.p_b.data, form.x_c.data,
                         form.T_env.data, form.eta.data)
    else:
        result = None

    return render_template('view.html', form=form, result=result)
if __name__ == '__main__':
    app.run(debug=False)

# http://127.0.0.1:5000/chaofan
