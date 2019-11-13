
class LayerConfig(object):

    _defaults = dict(c=3, a=0)

    def __init__(self, **kw):
        self._kw = LayerConfig._defaults.copy()
        
        self._kw.update(**kw)

    def __str__(self):
        return str(self._kw)

    def __getattr__(self, name):
        return self._kw.get(name)



x = LayerConfig(a=1, b=2)
y = LayerConfig(a=2, b=2, c=2)
z = LayerConfig()
print(x)
print(y)
print(z)


