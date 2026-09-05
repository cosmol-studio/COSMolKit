<script setup>
import { onBeforeUnmount, onMounted, ref, watch } from "vue";
import init, { Molecule } from "@cosmol-studio/cosmolkit";

const props = defineProps({ smiles: { type: String, default: "CCO" } });
const molecule = ref(null);
const error = ref("");
let ready = false;

function update() {
  if (!ready) return;
  try {
    const next = Molecule.fromSmiles(props.smiles);
    molecule.value?.free();
    molecule.value = next;
    error.value = "";
  } catch (value) {
    error.value = String(value);
  }
}

onMounted(async () => {
  try {
    await init();
    ready = true;
    update();
  } catch (value) {
    error.value = String(value);
  }
});

watch(() => props.smiles, update);
onBeforeUnmount(() => molecule.value?.free());
</script>

<template>
  <pre v-if="error">{{ error }}</pre>
  <output v-else-if="molecule">{{ molecule.toSmiles() }}</output>
  <p v-else>Loading COSMolKit...</p>
</template>
