def findCurve(prime, curveCofactor, twistCofactor, _A):
        Fp = GF(prime)
        A = _A
        while A < _A + 1000000000:
                # print(A)
                if(A -2) % 4 != 0:
                        A +=1
                        continue
                try:
                        E = EllipticCurve(Fp , [0, A, 0, 1, 0])
                except:
                        A +=1
                        continue

                groupOrder = E.order()
                if(groupOrder % curveCofactor != 0 or not is_prime(groupOrder // curveCofactor)):
                        A +=1
                        continue

                twistOrder = 2*(prime +1)-groupOrder
                if(twistOrder % twistCofactor != 0 or not is_prime(twistOrder // twistCofactor)):
                        A +=1
                        continue

                return E, A, 1, groupOrder, curveCofactor, groupOrder // curveCofactor

def find1Mod4(prime, curveCofactor, twistCofactor, A):
        assert((prime % 4) == 1)
        return findCurve(prime, curveCofactor, twistCofactor, A)

# Baby Jubjub in Montgomery form
###########################################################################
prime = 21888242871839275222246405745257275088548364400416034343698204186575808495617
Fp = GF(prime)
h = 8
A = 168688 # Start by 1 if you have time
EC, A, B, n, h, l = find1Mod4(prime, h, 4, A)
###########################################################################
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Step 1: Choice of Montgomery Equation")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("EC: ", EC)
print("Montgomery A: ", A)
print("Montgomery B: ", B)
print("Order of the curve: ", n)
print("Cofactor of the curve: ", h)
print("Largest prime divisor of the curves order: ", l)
###########################################################################

def findGenPoint(prime , A, EC, N):
    Fp = GF(prime)
    #counter = 1
    for uInt in range(1, 1000):
        #print(counter)
        #counter=counter + 1
        u = Fp(uInt)
        v2 = u^3 + A*u^2 + u
        if not v2.is_square():
            continue
        v = v2.sqrt()

        point = EC(u, v)
        pointOrder = point.order()
        if pointOrder == N:
            return point

def findBasePoint(EC, h, u, v):
    return h*EC(u, v)

# Generator and base points of Baby Jubjub in Montgomery form
gen_u, gen_v, gen_w = findGenPoint(prime, A, EC, n)
base_u, base_v, base_w = findBasePoint(EC, h, gen_u, gen_v)
##########################################################################
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Step 2: Choice of Generator and Base Points")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Montgomery Generator Point: ", gen_u, gen_v, gen_w)
print("Montgomery Base Point: ", base_u, base_v, base_w)
##########################################################################

def mont_to_ted (u, v, prime):
    Fp = GF(prime)
    x = Fp(u / v)
    y = Fp ((u -1)/(u +1))
    return (x, y)

def ted_to_mont(x, y, prime):
    Fp = GF(prime)
    u = Fp((1 + y)/(1 - y))
    v = Fp((1 + y )/((1 - y) * x))
    return (u, v)

def is_on_ted (x, y, prime, a, d):
    Fp = GF(prime)
    return Fp(a*(x **2) + y**2 - 1 - d*(x **2)*( y **2)) == 0

# Conversion of Baby Jubjub to twisted Edwards
a = Fp((A + 2) / B)
d = Fp((A - 2) / B)

# Check we have a safe twist and discriminant != 0
assert(not d.is_square())
assert(a*d*(a-d)!=0)

# Conversion of generator to twisted Edwards
gen_x, gen_y = mont_to_ted(gen_u, gen_v, prime)
assert (is_on_ted(gen_x, gen_y, prime, a, d))

# Sanity check : the inverse map returns the original point in Montgomery
u, v = ted_to_mont(gen_x, gen_y, prime)
assert (u == gen_u)
assert (v == gen_v)

# Conversion of base point to twisted Edwards
base_x, base_y = mont_to_ted(base_u, base_v, prime)
assert(is_on_ted(base_x, base_y, prime, a, d))

#########################################################
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Step 3: Transformation to Twisted Edwards")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Twisted Edwards a: ", a)
print("Twisted Edwards d: ", d)
print("Twisted Edwards Generator Point: ", gen_x, gen_y)
print("Twisted Edwards Base Point: ", base_x, base_y)
#########################################################

def scaling (a, d, prime):
    Fp = GF(prime)
    if Fp(-a).is_square():
        f = sqrt(Fp(-a));
        a_ = Fp(a/(f*f));
        d_ = Fp(d/(-a));
        if a_ == Fp(-1):
            a_ = -1
        else :
            a_ = a;
            d_ = a;
    return a_, d_, f

def ted_to_tedprime(x, y, prime, scaling_factor):
    Fp = GF(prime)
    x_ = Fp(x * (- scaling_factor))
    y_ = y;
    return(x_ , y_)

def tedprime_to_ted(x_, y_, prime, scaling_factor):
    Fp = GF(prime)
    x = Fp(x_ / (- scaling_factor))
    y = y_
    return(x, y)

def is_on_ted_prime(x, y, prime, a_, d_):
    Fp = GF( prime )
    return Fp(a_ *(x **2) + y **2 - 1 - d_ *(x **2)*( y **2)) == 0

# Conversion of E to Eopt
a_, d_, f = scaling(a, d, prime)

# Conversion of generator to Eopt
gen_x_prime, gen_y_prime = ted_to_tedprime(gen_x, gen_y, prime, f);
assert(is_on_ted_prime(gen_x_prime, gen_y_prime, prime, a_, d_))

# Sanity check : the inverse map returns the original point in E
u, v = tedprime_to_ted(gen_x_prime, gen_y_prime, prime, f)
assert(u == gen_x)
assert(v == gen_y)

# Conversion of base point to Eopt
base_x_prime , base_y_prime = ted_to_tedprime(base_x, base_y, prime, f);
assert(is_on_ted_prime(base_x_prime, base_y_prime, prime, a_ , d_))

# Sanity check : the inverse map returns the original point in E
u, v = tedprime_to_ted(base_x_prime, base_y_prime, prime, f)
assert(u == base_x)
assert(v == base_y)
#####################################################################
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Step 4: Optimization of Parameters")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Optimized Twisted Edwards a: ", a_)
print("Optimized Twisted Edwards d: ", d_)
print("Optimized Twisted Edwards scaling factor f: ", f)
print("Optimized Twisted Edwards Generator Point: ", gen_x_prime, gen_y_prime)
print("Optimized Twisted Edwards Base Point: ", base_x_prime, base_y_prime)
