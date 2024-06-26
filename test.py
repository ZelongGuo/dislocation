import dislocation as dl
import numpy as np 

mu = 3e10
nu = 0.25

#model = np.array([
#    440               ,3940             ,15,80,50,50, 45,    1,    1, 0, 
#	])

model = np.array([
    [440.58095043254673,3940.114839963042,15,80,50,50, 45, 0.01, 0.01, 0], 
	[440.58095043254673,3940.114839963042,15,80,50,50, 45, 0.01, 0.01, 0]])
# print(f"model: {model}")

obs = np.array([
    [454, 3943, 10],
    [454, 3943, 10],
    [454, 3943, 0],
    [454, 3943, 0],
    [454, 3943, 0],
    [454, 3943, 0],
    [454, 3943, 0],
    [454, 3943, 0]])
#obs = np.array([454, 3943, 0])
#print(f"obs: {obs}")

# length, width, depth, dip, strike, easting, northing, str-slip, dip-selip. opening


results = dl.okada_rect(obs, model, mu, nu)
#print(f"U = {results[0]}")
#print(f"D = {results[1]}")
print(f"S = {results[2]}")
print(f"E = {results[3]}")
print(f"flags = {results[4]}")



