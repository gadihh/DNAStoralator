from PyQt5.QtWidgets import QDoubleSpinBox


class SpinBoxCustom(QDoubleSpinBox):
    def __init__(self, *args, **kwargs):
        super(SpinBoxCustom, self).__init__(*args, **kwargs)

    def stepBy(self, steps: int) -> None:
        if steps == 1:
            self.setValue(self.value() * 10)
        if steps == -1:
            self.setValue(self.value() / 10)

    def textFromValue(self, v: float) -> str:
        no_zero = f'{v:g}'
        return str(no_zero)
