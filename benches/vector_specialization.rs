// benches/my_algo.rs
use criterion::{Criterion, criterion_group, criterion_main};
use std::{hint::black_box, marker::PhantomData, mem::ManuallyDrop, ops::Deref, ptr};

// From module `view`.
/// An immutable view to `T`.
pub struct View<'a, T> {
    inner: ManuallyDrop<T>,
    must_drop: bool,
    phantom: PhantomData<&'a T>,
}

impl<T> Drop for View<'_, T> {
    fn drop(&mut self) {
        if self.must_drop {
            unsafe { ManuallyDrop::drop(&mut self.inner) };
        }
    }
}

impl<T> Deref for View<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<'a, T> View<'a, T> {
    pub(crate) fn new(t: T, must_drop: bool) -> Self {
        Self {
            inner: ManuallyDrop::new(t),
            must_drop,
            phantom: PhantomData,
        }
    }
}

/// Vector type.
struct VecF64 {
    inner: *mut sys::gsl_vector,
}

impl VecF64 {
    fn new(n: usize) -> Self {
        let ptr = unsafe { sys::gsl_vector_calloc(n) };
        VecF64 { inner: ptr }
    }

    fn as_ptr(&self) -> *const sys::gsl_vector {
        self.inner
    }

    fn len(&self) -> usize {
        unsafe { (*self.as_ptr()).size }
    }

    fn stride(&self) -> usize {
        unsafe { (*self.as_ptr()).stride }
    }

    fn as_slice(&self) -> &[f64] {
        let data = unsafe { (*self.as_ptr()).data };
        unsafe { std::slice::from_raw_parts(data, self.len()) }
    }

    fn as_gsl_vector(x: &Self) -> View<'_, sys::gsl_vector> {
        let v = sys::gsl_vector {
            size: x.len(),
            stride: x.stride(),
            // `View` only offers read-only access so the cast is OK.
            data: VecF64::as_slice(x).as_ptr() as *mut _,
            block: ptr::null_mut(),
            owner: 0,
        };
        View::new(v, false)
    }

    fn as_gsl_vector_custom(x: &Self) -> View<'_, sys::gsl_vector> {
        let t = unsafe { *x.inner };
        View::new(t, false)
    }
}

fn bench_impls(c: &mut Criterion) {
    let v = VecF64::new(10);
    c.bench_function("ptr", |b| b.iter(|| black_box(&v).as_ptr()));
    c.bench_function("generic", |b| {
        b.iter(|| VecF64::as_gsl_vector(black_box(&v)))
    });
    c.bench_function("customized", |b| {
        b.iter(|| VecF64::as_gsl_vector_custom(black_box(&v)))
    });
}

criterion_group!(benches, bench_impls);
criterion_main!(benches);

// customized ≈ 2.35 ptr
// generic ≈ 1.25 customized
//
// However, we are speaking of a few hundred picoseconds!
