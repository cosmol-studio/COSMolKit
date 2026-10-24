//! Native scalar configuration projection. Unsupported shapes fail compilation.
//! No business-type map or duplicated native field/default declarations.
#[macro_export]
macro_rules! project_parameters {
    (@resume {$name:ident} @adapter $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($name, $native { $($(#[$attr])* $field: $ty),* });
    };
    (@resume {@flatten $name:ident, $native:path {$($parent:tt)*}, $field:ident} @$marker:ident $child:path { $($children:tt)* } methods {}) => {
        cosmolkit_tsify_macros::flatten_parameters!($name,$native{$($parent)*},$field,{$($children)*});
    };
    (@resume {@flatten $name:ident, $native:path {$($parent:tt)*} factory, $field:ident} @$marker:ident $child:path { $($children:tt)* } methods {}) => {
        cosmolkit_tsify_macros::flatten_parameters!($name,$native{$($parent)*} factory,$field,{$($children)*});
    };
    (@resume {$name:ident} @fields_factory $native:path { $($fields:tt)* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($name,$native{$($fields)*} factory);
    };
    (@resume {$name:ident} @native_factory $native:path { $($fields:tt)* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($name,$native{$($fields)*} factory);
    };
    (@resume {$name:ident} @native $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($name, $native { $($(#[$attr])* $field: $ty),* });
    };
    (@resume {$name:ident} @fields $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($name, $native { $($(#[$attr])* $field: $ty),* });
    };
}
