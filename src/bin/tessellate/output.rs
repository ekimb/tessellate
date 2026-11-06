use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::mpsc::{self, Receiver, SyncSender, TryRecvError, TrySendError};
use std::thread::{self, JoinHandle};
use tessellate::utils::create_with_capacity;

const SYNC_BYTES: u64 = 256 << 20;
const BACKLOG_BYTES: u64 = 1 << 30;

pub(super) struct OutputFile {
    inner: BufWriter<File>,
    requests: Option<SyncSender<u64>>,
    completed: Receiver<io::Result<u64>>,
    worker: Option<JoinHandle<()>>,
    pub bytes: u64,
    synced: u64,
    submitted: u64,
    next_checkpoint: u64,
    checkpoint_bytes: u64,
    backlog_bytes: u64,
}

fn stopped() -> io::Error {
    io::Error::other("output synchronization worker stopped")
}

impl OutputFile {
    pub fn new(path: &str, buffer_bytes: usize) -> Self {
        let inner = create_with_capacity(path, buffer_bytes);
        let file = inner
            .get_ref()
            .try_clone()
            .expect("cannot clone output file");
        Self::with_sync(inner, SYNC_BYTES, BACKLOG_BYTES, move || file.sync_data())
    }

    fn with_sync(
        inner: BufWriter<File>,
        checkpoint_bytes: u64,
        backlog_bytes: u64,
        mut sync: impl FnMut() -> io::Result<()> + Send + 'static,
    ) -> Self {
        assert!(checkpoint_bytes > 0 && checkpoint_bytes <= backlog_bytes);
        let (requests, pending) = mpsc::sync_channel(1);
        let (done, completed) = mpsc::channel();
        let worker = thread::spawn(move || {
            for position in pending {
                let result = sync().map(|()| position);
                let failed = result.is_err();
                if done.send(result).is_err() || failed {
                    break;
                }
            }
        });
        Self {
            inner,
            requests: Some(requests),
            completed,
            worker: Some(worker),
            bytes: 0,
            synced: 0,
            submitted: 0,
            next_checkpoint: checkpoint_bytes,
            checkpoint_bytes,
            backlog_bytes,
        }
    }

    fn poll(&mut self) -> io::Result<()> {
        loop {
            match self.completed.try_recv() {
                Ok(result) => self.synced = result?,
                Err(TryRecvError::Empty) => return Ok(()),
                Err(TryRecvError::Disconnected) => return Err(stopped()),
            }
        }
    }

    fn wait(&mut self) -> io::Result<()> {
        self.synced = self.completed.recv().map_err(|_| stopped())??;
        Ok(())
    }

    fn checkpoint(&mut self) -> io::Result<()> {
        self.inner.flush()?;
        self.poll()?;
        match self.requests.as_ref().unwrap().try_send(self.bytes) {
            Ok(()) => self.submitted = self.bytes,
            Err(TrySendError::Full(_)) => {}
            Err(TrySendError::Disconnected(_)) => return Err(stopped()),
        }
        self.next_checkpoint = self.bytes.saturating_add(self.checkpoint_bytes);
        Ok(())
    }
}

impl Write for OutputFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        if self.bytes - self.synced == self.backlog_bytes {
            self.poll()?;
            if self.bytes - self.synced == self.backlog_bytes {
                self.wait()?;
            }
        }
        let available = self.backlog_bytes - (self.bytes - self.synced);
        let n = self
            .inner
            .write(&buf[..buf.len().min(available as usize)])?;
        self.bytes += n as u64;
        if self.bytes >= self.next_checkpoint || self.bytes - self.synced == self.backlog_bytes {
            self.checkpoint()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()?;
        self.poll()?;
        if self.submitted < self.bytes {
            self.requests
                .as_ref()
                .unwrap()
                .send(self.bytes)
                .map_err(|_| stopped())?;
            self.submitted = self.bytes;
        }
        while self.synced < self.bytes {
            self.wait()?;
        }
        Ok(())
    }
}

impl Drop for OutputFile {
    fn drop(&mut self) {
        drop(self.requests.take());
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}
