__author__ = 'moore'

from flask import Flask
from views import *

app = Flask(__name__)
app.config.from_object('websiteconfig')


@app.errorhandler(404)
def not_found(error):
    return render_template('404.html'), 404


app.add_url_rule('/', view_func=index)

app.add_url_rule('/rc/', view_func=rc, methods=['GET'])
app.add_url_rule('/rcplot/', view_func=rcplot, methods=['POST'])
app.add_url_rule('/rcplot/<value>/', view_func=rcplot, methods=['POST'])

app.add_url_rule('/plot/<value1>/', view_func=plot, methods=['POST'])
app.add_url_rule('/plot/<value1>/<value2>/', view_func=plot, methods=['POST'])

app.add_url_rule('/lr/', view_func=lr, methods=['GET'])
app.add_url_rule('/lrplot/', view_func=lrplot, methods=['POST'])
app.add_url_rule('/lrplot/<value>/', view_func=lrplot, methods=['POST'])

app.add_url_rule('/lrc_series/', view_func=lrc_series, methods=['POST', 'GET'])
app.add_url_rule('/lrc_parallel/', view_func=lrc_parallel, methods=['POST', 'GET'])
app.add_url_rule('/lrc/<value>/', view_func=lrc, methods=['POST', 'GET'])
app.add_url_rule('/lrcplot/', view_func=lrcplot, methods=['POST', 'GET'])