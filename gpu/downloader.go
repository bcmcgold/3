package gpu

import (
	"github.com/barnex/cuda4/cu"
	"nimble-cube/core"
)

// Downloads data from GPU to host.
type Downloader struct {
	dev    RChan
	host   core.Chan
	bsize  int
	stream cu.Stream
}

func NewDownloader(devdata RChan, hostdata core.Chan) *Downloader {
	core.Assert(hostdata.Size() == devdata.Size())
	blocklen := core.Prod(core.BlockSize(hostdata.Size()))
	return &Downloader{devdata, hostdata, blocklen, 0} // TODO: block size
}

func (u *Downloader) Run() {
	core.Debug("run gpu.downloader with block size", u.bsize)
	LockCudaThread()
	u.stream = cu.StreamCreate()
	MemHostRegister(u.host.UnsafeData())

	for {
		in := u.dev.ReadNext(u.bsize)
		out := u.host.WriteNext(u.bsize)
		in.CopyDtoHAsync(out, u.stream)
		u.stream.Synchronize()
		u.dev.ReadDone()
		u.host.WriteDone()
	}
}
