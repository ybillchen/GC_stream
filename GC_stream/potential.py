# Licensed under BSD-3-Clause License - see LICENSE

class Potential(object):
	"""
	Time-dependent galaxy Potential
	"""

	def __init__(self, arg):
		super(Potential, self).__init__()
		self.arg = arg

	def generate_at(self, tlb):
		"""
		Interporlate potential at lookback time
		"""
		pass
		