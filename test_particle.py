import unittest
import numpy as np
from particle import Particle
from particle import read_particles

class TestParticle(unittest.TestCase):

    def setUp(self):
        self.test_particle = Particle("test",1,1,1,0,0,0,1,0)
    
    def test_particle_name(self):
        self.assertEqual(self.test_particle.name, "test")
    
    def test_particle_mass(self):
        self.assertEqual(self.test_particle.Mass, float(0))
    
    def test_particle_J(self):
        self.assertEqual(self.test_particle.J, 0)
    
    def test_particle_charm(self):
        self.assertEqual(self.test_particle.Charm, 0)

    def test_particle_strangeness(self):
        self.assertEqual(self.test_particle.Strange, 0)
    
    def test_particle_isospin(self):
        self.assertEqual(self.test_particle.Isospin, 1)
    
    def test_particle_isospin_charm(self):
        self.assertEqual(self.test_particle.Isospin_charm, 1)
    
    def test_read_particles(self):
        particles = read_particles("Particles/test.txt")
        self.assertEqual(len(particles), 2)




