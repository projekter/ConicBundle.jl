dll void cb_pscbundleparameters_destroy(PSCBundleParameters* self) {
  delete self;
}

dll int cb_pscbundleparameters_init(PSCBundleParameters* self, const BundleParameters* bp) {
  return self->init(*bp);
}

dll PSCBundleParameters* cb_pscbundleparameters_new() {
  return new PSCBundleParameters();
}

dll PSCBundleParameters* cb_pscbundleparameters_new2(const BundleParameters* bp) {
  return new PSCBundleParameters(*bp);
}

