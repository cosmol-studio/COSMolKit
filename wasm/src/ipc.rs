//! Borrowed distance matrices and numeric row overrides for graph information content.
use crate::{Molecule, errors::enum_error};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;

crate::enums::enum_value!(IpcMethod, "IpcMethod" {
    RdkitCompatible = "rdkit_compatible", ExactPolynomial = "exact_polynomial"
});

#[derive(serde::Deserialize, tsify::Tsify)]
#[serde(transparent)]
struct IpcRows(Vec<Vec<f64>>);

#[wasm_bindgen]
extern "C" {
    #[cfg(feature = "cap-matrices")]
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<number>> | DenseMatrix")]
    pub type IpcMatrixInput;
    #[cfg(not(feature = "cap-matrices"))]
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<number>>")]
    pub type IpcMatrixInput;

    #[cfg(feature = "cap-matrices")]
    #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
    fn borrow_matrix(
        visit: &mut dyn FnMut(&crate::matrix_calls::DenseMatrix),
        receiver: &JsValue,
        args: &js_sys::Array,
    ) -> Result<JsValue, JsValue>;
}

impl IpcMatrixInput {
    fn with_matrix<R>(
        &self,
        apply: impl FnOnce(ck::IpcDistanceMatrix<'_>) -> Result<R, JsValue>,
    ) -> Result<R, JsValue> {
        let raw: &JsValue = self.as_ref();
        if raw.is_instance_of::<js_sys::Array>() {
            let rows = crate::input::decode(tsify::Ts::<IpcRows>::new_unchecked(raw.clone()))?.0;
            return apply(rows.as_slice().into());
        }

        #[cfg(feature = "cap-matrices")]
        {
            crate::input::with_instance(
                raw,
                "DenseMatrix",
                |visit, receiver, args| borrow_matrix(&mut |matrix| visit(matrix), receiver, args),
                |matrix| apply((&matrix.0).into()),
            )
        }
        #[cfg(not(feature = "cap-matrices"))]
        {
            Err(js_sys::TypeError::new("expected an array of numeric rows").into())
        }
    }
}

#[wasm_bindgen]
impl Molecule {
    /// Total information content using an explicit polynomial algorithm.
    /// ExactPolynomial uses exact integer coefficients and floating-point entropy.
    #[wasm_bindgen(js_name = ipcWithMethod)]
    pub fn ipc_with_method(&self, method: IpcMethodInput) -> Result<f64, JsValue> {
        use crate::fields::Field;
        let method = ck::IpcMethod::input(method)?;
        self.0.ipc_with_method(method).map_err(IpcError::from_rust)
    }

    /// Average information content using an explicit polynomial algorithm.
    /// ExactPolynomial uses exact integer coefficients and floating-point entropy.
    #[wasm_bindgen(js_name = avgIpcWithMethod)]
    pub fn avg_ipc_with_method(&self, method: IpcMethodInput) -> Result<f64, JsValue> {
        use crate::fields::Field;
        let method = ck::IpcMethod::input(method)?;
        self.0
            .avg_ipc_with_method(method)
            .map_err(IpcError::from_rust)
    }

    pub fn ipc(&self, matrix: Option<IpcMatrixInput>) -> Result<f64, JsValue> {
        match matrix {
            Some(matrix) => self.ipc_with_distance_matrix(matrix),
            None => self.0.ipc().map_err(IpcError::from_rust),
        }
    }

    #[wasm_bindgen(js_name = ipcWithDistanceMatrix)]
    pub fn ipc_with_distance_matrix(&self, matrix: IpcMatrixInput) -> Result<f64, JsValue> {
        matrix.with_matrix(|m| {
            self.0
                .ipc_with_distance_matrix(m)
                .map_err(IpcError::from_rust)
        })
    }

    #[wasm_bindgen(js_name = avgIpc)]
    pub fn avg_ipc(&self, matrix: Option<IpcMatrixInput>) -> Result<f64, JsValue> {
        match matrix {
            Some(matrix) => self.avg_ipc_with_distance_matrix(matrix),
            None => self.0.avg_ipc().map_err(IpcError::from_rust),
        }
    }

    #[wasm_bindgen(js_name = avgIpcWithDistanceMatrix)]
    pub fn avg_ipc_with_distance_matrix(&self, matrix: IpcMatrixInput) -> Result<f64, JsValue> {
        matrix.with_matrix(|m| {
            self.0
                .avg_ipc_with_distance_matrix(m)
                .map_err(IpcError::from_rust)
        })
    }
}

enum_error!(IpcError, IpcErrorKind, ck::IpcError, {
    InvalidDistanceMatrixShape => ck::IpcError::InvalidDistanceMatrixShape { .. },
    InvalidBondEndpoint => ck::IpcError::InvalidBondEndpoint { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidDistanceMatrixShape { row, .. } => row, _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidDistanceMatrixShape { actual, .. } => Some(actual), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidDistanceMatrixShape { expected, .. } => Some(expected), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidBondEndpoint { bond, .. } => Some(bond), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn endpoint(&self) -> Option<String> { match *self.inner {
        ck::IpcError::InvalidBondEndpoint { endpoint, .. } => Some(endpoint.to_owned()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidBondEndpoint { atom, .. } => Some(atom), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner {
        ck::IpcError::InvalidBondEndpoint { atom_count, .. } => Some(atom_count), _ => None,
    } }
});
