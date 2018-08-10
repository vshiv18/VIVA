from varsome_api.client import VarSomeAPIClient, VarSomeAPIException
from varsome_api.models.variant import AnnotatedVariant
api_key = 'DqBQzz@0ON&?y%UO6ECrKxAj4QP4SDU7k@hROW&5'
api = VarSomeAPIClient(api_key)
revel = {}
def single_lookup(id):
	#a single lookup in the Varsome API, returns a single annotated variant wrapper
	return api.lookup(id, params={'add-all-data':1})
def batch_lookup(lst):
	#given a file with variant ids on each line
	#returns a list of annotated variant wrappers
	return api.batch_lookup(lst, params={'add-all-data':1})
def revel_score(variant):
	#given an annotated variant object
	#returns the revel score
	return revel[variant.pos]

