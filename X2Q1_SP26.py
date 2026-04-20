# region imports
import numpy as np
import math
from scipy.optimize import fsolve
import random as rnd
# endregion

# Please note that Perplexity was used to help with debugging and writing

# region class definitions
class UC():  # a units conversion class
    def __init__(self):
        """
        This unit converter class is useful for the pipe network and perhaps other problems.
        The strategy is (number in current units)*(conversion factor)=(number desired units), for instance:
            1(ft)*(self.ft_to_m) = 1/3.28084 (m)
            1(in^2)*(self.in2_to_m2) = 1*(1/(12*3.28084))**2 (m^2)
        """

    #region class constants
    #we have not used these in class yet, but I think they are not too mysterious.
    ft_to_m = 1 / 3.28084
    ft2_to_m2 = ft_to_m ** 2
    ft3_to_m3 = ft_to_m ** 3
    ft3_to_L = ft3_to_m3 * 1000
    L_to_ft3 = 1 / ft3_to_L
    in_to_m = ft_to_m / 12
    m_to_in = 1 / in_to_m
    in2_to_m2 = in_to_m ** 2
    m2_to_in2 = 1 / in2_to_m2
    g_SI = 9.80665  # m/s^2
    g_EN = 32.174  # 32.174 ft/s^2
    gc_EN = 32.174  # lbm*ft/lbf*s^2
    gc_SI = 1.0  # kg*m/N*s^2
    lbf_to_kg = 1 / 2.20462
    lbf_to_N = lbf_to_kg * g_SI
    pa_to_psi = (1 / (lbf_to_N)) * in2_to_m2
    #endregion

    @classmethod  # this notation allows this method to be directly used from the class by UC.viscosityEnglishToSI
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Converts between lb*s/ft^2 and Pa*s
        :param mu: the viscosity in english units
        :param toSI:  True assumes english in, False assumes SI in
        :return: the viscosity in Pa*s if toSI=True, lb*s/ft^2 if toSI=False
        """
        #(lb*s)/ft^2*((3.3 ft/m)^2)*(1kg/2.2lb)*(9.81m/s^2)->(Pa*s)
        cf=(1/cls.ft2_to_m2)*(cls.lbf_to_kg)*cls.g_SI
        return mu*cf if toSI else mu/cf

    @classmethod
    def densityEnglishToSI(cls, rho, toSI=True ):
        """
        Converts between lb/ft^3 and kg/m^3
        :param rho: specific weight or density
        :param toSI:  True assumes english in, False assumes SI in
        :return: density in SI or EN
        """
        #(lb/ft^3)*((3.3ft/m)^3)*(1kg/2.2lb) -> kg/m^3
        cf=cls.lbf_to_kg/cls.ft3_to_m3
        return rho*cf if toSI else rho/cf

    @classmethod
    def head_to_pressure(cls, h, rho, SI=True):
        """
        Convert from height of column of fluid to pressure in consistent units
        :param h: head in height of fluid (in or m)
        :return: pressure in (psi or Pa)
        """
        if SI: # p = rho*g*h = g*cf
            cf=rho*cls.g_SI/cls.gc_SI  # kg*m/m^3*s^2
            return h*cf
        else:  # p = rho*g*h = g*cf (h in in)
            cf = rho*cls.g_EN/cls.gc_EN*(1/12)**2  # (lbm*ft/ft^3*s^2)(lbf*s^2/lbm*ft)(ft^2/in^2)
            return h*cf
        #convert m of water to psi
        #(m)*(3.3*12in/m)*rho(kg/m^3)*(2.2lb/kg)*(1m/(3.3*12in))^3
        psi=p*cls.rho*2.2/((3.3*12)**2)
        return psi

    @classmethod
    def m_to_psi(cls, h, rho):
        """
        For converting from height of fluid to psi
        :param h: height of fluid in m
        :param rho: density of fluid in kg/m^3
        :return: pressure in psi
        """
        return cls.head_to_pressure(h,rho)*cls.pa_to_psi
    
    @classmethod
    def psi_to_m(cls, p, rho):
        """
        For converting from psi to height of fluid.
        first convert psi to pa
        :param p: pressure in psi
        :param rho: density of fluid in kg/m^3
        :return: height of fluid in m
        """
        pa=p/cls.pa_to_psi
        h= pa/(rho*cls.g_SI)
        return h

class Fluid():
    def __init__(self, mu=0.00089, rho=1000, SI=True):
        '''
        default properties are for water in the SI system
        :param mu: dynamic viscosity in Pa*s -> (kg*m/s^2)*(s/m^2) -> kg/(m*s) or (lb*s/ft^2)
        :param rho: density in kg/m^3 or (lb/ft^3)
        :param SI: tells constructor if unit conversion is needed from english to si if SI==False
        '''
        self.mu= mu if SI==True else UC.viscosityEnglishToSI(mu)
        self.rho= rho if SI==True else UC.densityEnglishToSI(rho)
        self.nu= self.mu/self.rho # calculate the kinematic viscosity in units of m^2/s

class Node():
    def __init__(self, Name='a', Pipes=[], ExtFlow=0):
        '''
        A node in a pipe network.
        :param Name: name of the node
        :param Pipes: a list/array of pipes objects connected to this node
        :param ExtFlow: any external flow into (+) or out (-) of this node in L/s
        '''
        self.name = Name
        self.pipes = Pipes
        self.extFlow = ExtFlow
        self.QNet=0
        self.P = 0 # the pressure head at the node in m of fluid
        self.oCalculated=False

    def getNetFlowRate(self):
        '''
        Calculates the net flow rate into this node in L/s
        # :return: the net flow rate into this node
        '''
        Qtot = self.extFlow  #count the external flow first
        for p in self.pipes:
            #retrieves the pipe flow rate (+) if into node (-) if out of node.  see class for pipe.
            Qtot += p.getFlowIntoNode(self.name)
        self.QNet=Qtot
        return self.QNet

    def setExtFlow(self, E, SI=True):
        """
        Sets the external flow rate for the node.  SI=False causes a conversion from ft^3/s to L/s
        :param E: External (volumetric) flow rate
        :param SI: boolean, False means units for E are in ft^3/s.  True means L/s
        :return: nothing
        """
        self.extFlow=E if SI else E*UC.ft3_to_L

class Loop():
    def __init__(self, Name='A', Pipes=[]):
        '''
        Defines a loop in a pipe network.  Note: the pipes must be listed in order.  The traversal of a pipe loop
        will begin at the start node of Pipe[0] and move in the positive direction of that pipe.  Hence, loops
        can be either CW or CCW traversed, depending on which pipe you start with.  Should work fine either way.
        :param Name: name of the loop
        :param Pipes: a list/array of pipes in this loop
        '''
        self.name=Name
        self.pipes=Pipes

    def getLoopHeadLoss(self):
        '''
        Calculates the net head loss as I traverse around the loop, in m of fluid.
        :return:
        '''
        deltaP=0 #initialize to zero
        startNode=self.pipes[0].startNode #begin at the start node of the first pipe
        for p in self.pipes:
            # calculates the head loss in the pipe considering loop traversal and flow directions
            phl=p.getFlowHeadLoss(startNode)
            deltaP+=phl
            startNode=p.endNode if startNode!=p.endNode else p.startNode #move to the next node

        return deltaP

class Pipe():
    def __init__(self, Start='A', End='B',L=100, D=200, r=0.00025, fluid=Fluid(), SI=True):
        '''
        Defines a generic pipe with orientation from lowest letter to highest, alphabetically.
        :param Start: the start node (string)
        :param End: the end node (string)
        :param L: the pipe length in m (float)
        :param D: the pipe diameter in mm (float)
        :param r: the pipe roughness in m  (float)
        :param fluid:  a Fluid object (typically water)
        :param SI: if SI==False, need to convert len, roughness from ft to m and dia from in to m
        '''
        # from arguments given in constructor
        
        self.startNode = min(Start.lower(), End.lower()) #force the alphabetically lower node name to be the pipe's start node
        self.endNode = max(Start.lower(),End.lower()) #makes sure to use the highest letter for the endNode
        self.length = L if SI else UC.ft_to_m*L
        self.rough = r if SI else UC.ft_to_m*r
        self.fluid = fluid #the fluid in the pipe

        # other calculated properties
        self.d=D/1000.0 if SI else UC.in_to_m*D #diameter in m
        self.relrough = self.rough / self.d #compute pipe cross-sectional area in m^2
        self.A = math.pi * self.d ** 2 / 4.0  #calculate pipe cross sectional area for easy use later
        self.Q = 10 #working in units of L/s, just an initial guess
        self.vel = self.V()  #calculate the initial velocity of the fluid
        self.reynolds = self.Re() #calculate the initial reynolds number
        self.hl = 0  #store the frictional head loss when calculated (m of fluid)

    def V(self):
        '''
        Calculate average velocity in the pipe for volumetric flow self.Q
        :return:the average velocity in m/s
        '''
        self.vel = (self.Q / 1000.0) / self.A # the average velocity is Q/A (be mindful of units)
        return self.vel

    def Re(self):
        '''
        Calculate the reynolds number under current conditions.
        :return:
        '''
        self.reynolds = self.fluid.rho * self.V() * self.d / self.fluid.mu # Re=rho*V*d/mu, be sure to use V() so velocity is updated.
        return self.reynolds

    def FrictionFactor(self):
        """
        This function calculates the friction factor for a pipe based on the
        notion of laminar, turbulent and transitional flow.
        :return: the (Darcy) friction factor
        """
        # update the Reynolds number and make a local variable Re
        Re = self.Re()
        rr = self.relrough
        # to be used for turbulent flow
        # laminar if truly laminar
        if Re < 2000:
            ff = 64.0 / Re
        else:
            # turbulent – Jain explicit form for Colebrook-White
            # f = 0.25 / [log10( (e/3.7D) + 5.74 / Re^0.9 )]^2
            ff = 0.25 / (math.log10(rr / 3.7 + 5.74 / (Re ** 0.9)) ** 2)

        return ff

    def frictionHeadLoss(self):  # calculate headloss through a section of pipe in m of fluid
        '''
        Use the Darcy-Weisbach equation to find the head loss through a section of pipe.
        DeltaP=f*(L/d)*(rho*V^2)/2
        Note:  the headloss should always be a positive number.
        '''
        g = 9.81  # m/s^2
        ff = self.FrictionFactor()
        v = self.V()  # recompute velocity from current Q
        self.hl = ff * (self.length / self.d) * (v ** 2) / (2 * g) # calculate the head loss in m of water
        return self.hl

    def getFlowHeadLoss(self, s):
        '''
        Calculate the head loss for the pipe.
        :param s: the node i'm starting with in a traversal of the pipe
        :return: the signed headloss through the pipe in m of fluid
        '''
        #while traversing a loop, if s = startNode I'm traversing in same direction as positive pipe
        nTraverse= 1 if s==self.startNode else -1
        #if flow is positive sense, scalar =1 else =-1
        nFlow=1 if self.Q >= 0 else -1
        return nTraverse*nFlow*self.frictionHeadLoss()

    def Name(self):
        '''
        Gets the pipe name.
        :return: pipe name (e.g., 'a-b')
        '''
        return self.startNode+'-'+self.endNode

    def oContainsNode(self, node):
        #does the pipe connect to the node?
        return self.startNode == node or self.endNode == node

    def printPipeFlowRate(self, SI=True):
        
        q_units = 'L/s' if SI else 'cfs'
        q=self.Q if SI else self.Q*UC.L_to_ft3
        print('The flow in segment {} is {:0.2f} ({}) and Re={:.1f}'.format(self.Name(),q,  q_units,self.reynolds))

    def printPipeHeadLoss(self, SI=True):
        """
        Print properties of pipe and head loss
        :param SI: if True, prints in mm of fluid, if False prints inches of fluid
        :return:
        """
        cfd=1000 if SI else UC.m_to_in  # conversion factor for diameter
        unitsd='mm' if SI else 'in'  # units for diameter
        cfL=1 if SI else 1/UC.ft_to_m # conversion factor for length
        unitsL='m' if SI else 'in'  # units for length
        cfh=cfd
        units_h=unitsd
        print("head loss in pipe {} (L={:.2f} {}, d={:.2f} {}) is {:.2f} {} of water".format(self.Name(),self.length*cfL,unitsL, self.d*cfd, unitsd,self.hl*cfh, units_h))

    def getFlowIntoNode(self, n):
        '''
        determines the flow rate into node n
        :param n: a node object
        :return: +/-Q
        '''
        if n == self.startNode:
            return -self.Q
        return self.Q

class PipeNetwork():
    def __init__(self, Pipes=[], Loops=[], Nodes=[], fluid=Fluid()):
        '''
        The pipe network is built from pipe, node, loop, and fluid objects.
        :param Pipes: a list of pipe objects
        :param Loops: a list of loop objects
        :param Nodes: a list of node objects
        :param fluid: a fluid object
        '''
        self.loops=Loops
        self.nodes=Nodes
        self.Fluid=fluid
        self.pipes=Pipes

    def findFlowRates(self):
        '''
        a method to analyze the pipe network and find the flow rates in each pipe
        given the constraints of: i) no net flow into a node and ii) no net pressure drops in the loops.
        :return: a list of flow rates in the pipes
        '''
        #see how many nodes and loops there are, this is how many equation results I will return
        N = len(self.pipes)
        # build an initial guess for flow rates in the pipes.
        # note that I only have 10 pipes, but need 11 variables because of the degenerate node equation at b.
        Q0 = np.full(N, 10.0)
        def fn(q):
            """
            This is used as a callback for fsolve.  The mass continuity equations at the nodes and the loop equations
            are functions of the flow rates in the pipes.  Hence, fsolve will search for the roots of these equations
            by varying the flow rates in each pipe.
            :param q: an array of flowrates in the pipes + 1 extra value b/c of node b
            :return: L an array containing flow rates at the nodes and  pressure losses for the loops
            """
            #Set all the nodes to 0 pressure and uncalculated
            for n in self.nodes:
                n.P = 0.0
                n.oCalculated = False
            #update the flow rate in each pipe object
            for i in range(N):
                self.pipes[i].Q = q[i]  # set volumetric flow rate from input argument q

            # node mass-balance equations
            qNet = self.getNodeFlowRates()  # length 10
            # drop one redundant node equation, e.g. node 'b'
            # nodes are built in order a, h, b, c, d, e, f, i, g, j in your buildNodes()
            # so find index for 'b' and remove it
            node_names = [n.name for n in self.nodes]
            idx_b = node_names.index('b')
            node_eqs = qNet[:idx_b] + qNet[idx_b + 1:]

            # loop head-loss equations (length 4)
            loop_eqs = self.getLoopHeadLosses()

            # total: 9 node eqs + 4 loop eqs = 13
            return node_eqs + loop_eqs

        #using fsolve to find the flow rates
        FR = fsolve(fn, Q0)

        # write solution back to pipes
        for i in range(len(self.pipes)):
            self.pipes[i].Q = FR[i]

        # recompute node flows and loop losses for reporting
        self.getNodeFlowRates()
        self.getLoopHeadLosses()

        return FR

    def getNodeFlowRates(self):
        #each node object is responsible for calculating its own net flow rate
        qNet = [n.getNetFlowRate() for n in self.nodes]
        return qNet

    def getLoopHeadLosses(self):
        #each loop object is responsible for calculating its own net head loss
        lhl=[l.getLoopHeadLoss() for l in self.loops]
        return lhl

    def getNodePressures(self, knownNodeP, knownNode):
        '''
        Calculates the pressures at the nodes by traversing the loops.  For this to work,
        I must traverse the nodes in the proper order, so that the start node of each loop except
        for the first one has been calculated before traversing the loop.
        :return:
        '''
        #step 1:  set all the node pressures to 0 and mark as not-calculated
        for n in self.nodes:
            n.P=0.0
            n.oCalculated=False

        #step 2:  visit nodes by traversing loops and set pressure of node at end of each pipe if not already calculated
        for l in self.loops:
            startNode = l.pipes[0].startNode  # begin at the start node of the first pipe
            n = self.getNode(startNode)  # get the startNode object
            CurrentP = n.P
            n.oCalculated=True
            for p in l.pipes:
                # calculates the head loss in the pipe considering loop traversal and flow directions
                phl = p.getFlowHeadLoss(startNode)
                CurrentP -= phl #subtract headloss to get positive pressure (head)
                startNode = p.endNode if startNode != p.endNode else p.startNode  # move to the next node
                n=self.getNode(startNode)  # get the node object
                n.P=CurrentP
                # n.oCalculated=True
        #step 3:  each node now has a P in  headloss, now set the knownNode pressure
        kn=self.getNode(knownNode)
        # presumably I started at some other node as my datum, so the known node should have a
        # headloss to all the other nodes
        deltaP = knownNodeP-1.0*kn.P  # will set head loss at known node to zero
        for n in self.nodes:
            n.P = n.P+deltaP  #sets positive pressure at nodes

    def getPipe(self, name):
        #returns a pipe object by its name
        for p in self.pipes:
            if name == p.Name():
                return p

    def getNodePipes(self, node):
        #returns a list of pipe objects that are connected to the node object
        l=[]
        for p in self.pipes:
            if p.oContainsNode(node):
                l.append(p)
        return l

    def nodeBuilt(self, node):
        #determines if I have already constructed this node object (by name)
        for n in self.nodes:
            if n.name==node:
                return True
        return False

    def getNode(self, name):
        #returns one of the node objects by name
        for n in self.nodes:
            if n.name==name:
                return n

    def buildNodes(self):
        #automatically create the node objects by looking at the pipe ends
        for p in self.pipes:
            if self.nodeBuilt(p.startNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.startNode,self.getNodePipes(p.startNode)))
            if self.nodeBuilt(p.endNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.endNode,self.getNodePipes(p.endNode)))

    def printPipeFlowRates(self, SI=True):
        for p in self.pipes:
            p.printPipeFlowRate(SI=SI)

    def printNetNodeFlows(self, SI=True):
        for n in self.nodes:
            Q=n.QNet if SI else n.QNet*UC.L_to_ft3
            units='L/S' if SI else 'cfs'
            print('net flow into node {} is {:0.2f} ({})'.format(n.name, Q, units))

    def printLoopHeadLoss(self, SI=True):
        cf=UC.m_to_psi(1,self.pipes[0].fluid.rho)
        units='m of water' if SI else 'psi'
        for l in self.loops:
            hl=l.getLoopHeadLoss()
            hl = hl if SI else hl*cf
            print('head loss for loop {} is {:0.2f} ({})'.format(l.name, hl, units))

    def printPipeHeadLoss(self, SI=True):
        cf=UC.m_to_in
        for p in self.pipes:
            p.printPipeHeadLoss(SI=SI)
            # hl=p.HeadLoss()
            # hl=hl if SI else hl*cf
            # l=p.length if SI else p.length/UC.ft_to_m
            # d=p.diam if SI else p.diam/UC.in_to_m
            # p_units='m of water' if SI else 'in of water'
            # print('head loss in pipe {} (L={:0.2f}, d={:0.3f}) is {:0.2f} {}'.format(p.Name(), l, d, hl, p_units))

    def printNodePressures(self, SI=True):
        pUnits='m of water' if SI else 'psi'
        cf = 1.0 if SI else UC.m_to_psi(1,self.Fluid.rho)
        for n in self.nodes:
            p = n.P*cf
            print('Pressure at node {} = {:0.2f} {}'.format(n.name,p,pUnits))
# endregion

# region function definitions
def main():
    '''
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc. (see problem statement)
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    '''
    #instantiate a Fluid object to define the working fluid as water
    #water= Fluid(mu=0.00089, rho=1000.0)
    SIUnits = False
    water = Fluid(mu=20.50e-6, rho=62.3, SI=SIUnits) # instantiate a fluid object

    # note:  these are surface roughness in ft.  We must calculate relative roughness for each pipe later
    r_CI=0.00085 #ft roughness for cast iron
    r_CN=0.003 #ft roughness for concrete

    #instantiate a new PipeNetwork object
    PN=PipeNetwork()
    PN.Fluid=water
    #add Pipe objects to the pipe network (see constructor for Pipe class)
    PN.pipes.append(Pipe('a','b',1000, 18, r_CN, water, SI=SIUnits))
    # Add all the necessary pipes to the network
    PN.pipes.append(Pipe('a', 'h', 1600, 24, r_CN, water, SI=SIUnits))  # left vertical pipe from a to h
    PN.pipes.append(Pipe('b', 'c', 500, 18, r_CN, water, SI=SIUnits))  # top horizontal pipe from b to c
    PN.pipes.append(Pipe('b', 'e', 800, 16, r_CI, water, SI=SIUnits))  # upper middle vertical pipe from b to e
    PN.pipes.append(Pipe('c', 'd', 500, 18, r_CN, water, SI=SIUnits))  # top horizontal pipe from c to d
    PN.pipes.append(Pipe('c', 'f', 800, 16, r_CI, water, SI=SIUnits))  # upper right vertical pipe from c to f
    PN.pipes.append(Pipe('d', 'g', 800, 16, r_CI, water, SI=SIUnits))  # far right vertical pipe from d to g
    PN.pipes.append(Pipe('e', 'f', 500, 12, r_CI, water, SI=SIUnits))  # middle horizontal pipe from e to f
    PN.pipes.append(Pipe('h', 'i', 1000, 24, r_CN, water, SI=SIUnits))  # bottom horizontal pipe from h to i
    PN.pipes.append(Pipe('f', 'g', 500, 12, r_CI, water, SI=SIUnits))  # middle horizontal pipe from f to g
    PN.pipes.append(Pipe('e', 'i', 800, 18, r_CN, water, SI=SIUnits))  # lower middle vertical pipe from e to i
    PN.pipes.append(Pipe('g', 'j', 800, 18, r_CN, water, SI=SIUnits))  # lower right vertical pipe from g to j
    PN.pipes.append(Pipe('i', 'j', 1000, 24, r_CN, water, SI=SIUnits))  # bottom horizontal pipe from i to j

    #add Node objects to the pipe network by calling buildNodes method of PN object
    PN.buildNodes()

    #update the external flow of certain nodes
    #update the external flow of inlet node
    PN.getNode('h').setExtFlow(10,SI=SIUnits)
    #set external flows at required nodes
    PN.getNode('e').setExtFlow(-3, SI=SIUnits)
    PN.getNode('f').setExtFlow(-5, SI=SIUnits)
    PN.getNode('d').setExtFlow(-2, SI=SIUnits)

    #add Loop objects to the pipe network
    # A: big left loop
    PN.loops.append(Loop('A', [
        PN.getPipe('a-b'),
        PN.getPipe('b-e'),
        PN.getPipe('e-i'),
        PN.getPipe('h-i'),
        PN.getPipe('a-h')
    ]))
    # add the remaining loops in the network
    # B: upper middle rectangle
    PN.loops.append(Loop('B', [
        PN.getPipe('b-c'),
        PN.getPipe('c-f'),
        PN.getPipe('e-f'),
        PN.getPipe('b-e')
    ]))

    # C: upper right rectangle
    PN.loops.append(Loop('C', [
        PN.getPipe('c-d'),
        PN.getPipe('d-g'),
        PN.getPipe('f-g'),
        PN.getPipe('c-f')
    ]))

    # D: lower central loop (note: e–f–i)
    PN.loops.append(Loop('D', [
        PN.getPipe('e-f'),
        PN.getPipe('f-g'),
        PN.getPipe('g-j'),
        PN.getPipe('i-j'),
        PN.getPipe('e-i')
    ]))

    #call the findFlowRates method of the PN (a PipeNetwork object)
    PN.findFlowRates()
    knownP = UC.psi_to_m(80, water.rho)
    PN.getNodePressures(knownNode='h', knownNodeP=knownP)

    #get output
    PN.printPipeFlowRates(SI=SIUnits)
    print()
    print('Check node flows:')
    PN.printNetNodeFlows(SI=SIUnits)
    print()
    print('Check loop head loss:')
    PN.printLoopHeadLoss(SI=SIUnits)
    print()
    PN.printPipeHeadLoss(SI=SIUnits)
    print()
    PN.printNodePressures(SI=SIUnits)
# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion

