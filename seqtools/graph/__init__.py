import uuid

class Node:
   """one point on the graph contains 'things'"""
   def __init__(self,payload=None,payload_list=[]):
      self._payload = []
      if payload is not None: self._payload.append(payload)
      if len(payload_list) > 0: self._payload = payload_list
      self.id = str(uuid.uuid4())
   @property
   def payload(self):
      if len(self._payload) == 0: return None 
      return self._payload[0]
   @property
   def payload_list(self): return self._payload
class Edge:
   """one edge on the graph contains 'things'"""
   def __init__(self,n1,n2,payload=None,payload_list=[]):
      self.node1 = n1
      self.node2 = n2
      self._payload = []
      if payload is not None: self._payload.append(payload)
      if len(payload_list) > 0: self._payload = payload_list
      self.id = str(uuid.uuid4())
   @property
   def payload(self):
      if len(self._payload) == 0: return None 
      return self._payload[0]
   @property
   def payload_list(self): return self._payload

class Graph(object):
   """Generic directed graph"""
   def __init__(self):
      self._nodes = {}
      self._edges = {}
      self._p2c = {} #parent to child by node id
      self._c2p = {} #child to parent by node id
   def __str__(self):
      o = ''
      o += str(len(self.roots))+" roots\n"
      o += str(len(self._nodes))+" nodes\n"
      o += str(len(self._edges))+" edges\n"
      o += str(sum([len(x.payload_list) for x in self._nodes.values()]))+" payload items\n"
      return o
   @property
   def nodes(self):  return self._nodes.values()

   @property
   def edges(self):  return self._edges.values()

   @property
   def roots(self):
      """get the nodes with no children"""
      return [x for x in self._nodes.values() if x.id not in self._c2p]

   def get_root_graph(self,root):
      """Return back a graph containing just the root and children"""
      children = self.get_children(root)
      g = Graph()
      nodes = [root]+children
      for node in nodes: g.add_node(node)
      node_ids = [x.id for x in nodes]
      edges = [x for x in self._edges.values() if x.node1.id in node_ids and x.node2.id in node_ids]
      for e in edges: g.add_edge(e)
      return g

   def get_children(self,root,_visited=[]):
      _visited = _visited[:]
      if root.id not in self._p2c: return _visited
      children = [self._nodes[x] for x in self._p2c[root.id].keys()]
      for c in children:
         _visited.append(c)
         _visited = self.get_children(c,_visited)
      return _visited

   def add_node(self,node):
      if node.id in self._nodes: return 
      self._nodes[node.id] = node

   def add_edge(self,edge):
      if edge.id in self._edges: return
      self.add_node(edge.node1)
      self.add_node(edge.node2)
      if edge.node1.id not in self._p2c:
         self._p2c[edge.node1.id] = {}
      if edge.node2.id not in self._p2c[edge.node1.id]:
         self._p2c[edge.node1.id][edge.node2.id] = {}
      self._p2c[edge.node1.id][edge.node2.id][edge.id] = edge

      if edge.node2.id not in self._c2p:
         self._c2p[edge.node2.id] = {}
      if edge.node1.id not in self._c2p[edge.node2.id]:
         self._c2p[edge.node2.id][edge.node1.id] = {}
      self._c2p[edge.node2.id][edge.node1.id][edge.id] = edge

      self._edges[edge.id] = edge
   def merge_cycles(self):
      """Work on this graph and remove cycles, with nodes containing concatonated lists of payloads"""
      while True:
         ### remove any self edges
         own_edges = self.get_self_edges()
         if len(own_edges) > 0:
            for e in own_edges: self.remove_edge(e)
         c = self.find_cycle()
         if not c: return
         keep = c[0]
         remove_list = c[1:]
         for n in remove_list: self.move_edges(n,keep)
         for n in remove_list: keep.payload_list += n.payload_list
         for n in remove_list: self.remove_node(n)

   def get_self_edges(self):
      edges = []
      for e in self._edges.values():
         if e.node1.id == e.node2.id: edges.append(e)
      return edges

   def remove_node(self,node):
      """remove the node"""
      if node.id not in self._nodes: return      
      """find edges to remove"""
      edges = set()
      for e in self._edges.values():
         if e.node1.id == node.id: edges.add(e.id)
         if e.node2.id == node.id: edges.add(e.id)
      edges = [self._edges[x] for x in list(edges)]
      for e in edges: self.remove_edge(e)
      del self._nodes[node.id]

   def remove_edge(self,edge):
      """Remove the edge"""
      if edge.id not in self._edges: return # its not in the graph
      del self._p2c[edge.node1.id][edge.node2.id][edge.id]
      if len(self._p2c[edge.node1.id][edge.node2.id].keys()) == 0:
         del self._p2c[edge.node1.id][edge.node2.id]
      if len(self._p2c[edge.node1.id].keys()) == 0:
         del self._p2c[edge.node1.id]

      del self._c2p[edge.node2.id][edge.node1.id][edge.id]
      if len(self._c2p[edge.node2.id][edge.node1.id].keys()) == 0:
         del self._c2p[edge.node2.id][edge.node1.id]
      if len(self._c2p[edge.node2.id].keys()) == 0:
         del self._c2p[edge.node2.id]

      del self._edges[edge.id]

   def move_edges(self,n1,n2):
      """Move edges from node 1 to node 2
 
         Not self edges though

         Overwrites edges
      """
      #Traverse edges to find incoming with n1
      incoming = []
      for e in self._edges.values():
         if e.node2.id == n1.id: incoming.append(e)
      #Traverse edges to find outgoing from n1
      outgoing = []
      for e in self._edges.values():
         if e.node1.id == n1.id: outgoing.append(e)
      #Make new edges to the new target
      for e in incoming:
         if e.node1.id == n2.id: continue # skip self
         newedge = Edge(e.node1,n2,payload_list=n2.payload_list+n1.payload_list)
         self.add_edge(newedge)
      for e in outgoing:
         if e.node2.id == n2.id: continue # skip self
         newedge = Edge(n2,e.node2,payload_list=n2.payload_list+n1.payload_list)
         self.add_edge(newedge)
      #now remove the edges that got transfered
      for e in incoming: self.remove_edge(e)
      for e in outgoing: self.remove_edge(e)

   def find_cycle(self):
      """greedy search for a cycle"""
      for node in self.nodes:
         cyc = self._follow_children(node)
         if len(cyc) > 0:
            return [self._nodes[x] for x in cyc]
      return None

   def _follow_children(self,node,parents=[]):
      parents = parents[:]
      if node.id in parents:
         while parents[0] != node.id: parents.pop(0)
         return parents
      parents.append(node.id)
      if node.id not in self._p2c: 
         return []
      children = [self._nodes[x] for x in self._p2c[node.id].keys()]
      for c in children:
         v = self._follow_children(c,parents[:])
         if len(v) > 0: return v
      return []
