import hashlib

__all__ = ["HashFile"]

class HashFile:
    def __init__(self, hasher="md5"):
        self.handler = None
        self.ashexstr = True
        self.blocksize = 16777216 # 2**24 = 16Mo
        self.internally_opened = False
        try:
            self.hasher = hashlib.__getattribute__(str(hasher))()
        except AttributeError as e:
            print(f"{e}")
            del self

    def _file_as_block(self):
        with self.handler:
            block = self.handler.read(self.blocksize)
            while len(block) > 0:
                yield block
                block = self.handler.read(self.blocksize)
    
    def hash_file(self, handler, blocksize=None):
        if isinstance(handler, str):
            self.handler = open(handler, "rb")
            self.internally_opened = True
        else:
            self.handler = handler

        if blocksize:
            self.blocksize = blocksize

        if self.handler.mode != "rb":
            raise ValueError("File must be open in read bytes mode !")

        for block in self._file_as_block():
            self.hasher.update(block)

        if self.internally_opened:
            self.handler.close()
            self.internally_opened = False

        return self.hasher.hexdigest() if self.ashexstr else self.hasher.digest()
    
