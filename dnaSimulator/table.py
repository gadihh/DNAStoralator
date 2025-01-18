from PyQt5.QtGui import *
from worker import *
from index_clustering import Clustering


class CellData:
    """
    Wrapper for data that is being stored in a cell. Some info is for maintenance and not for display,
    which will be stored in a hidden field.
    """

    def __init__(self, data, hidden=None):
        self.data = data
        if hidden is not None:
            self.hidden = hidden


class FilterModel(QSortFilterProxyModel):
    """
    Filtering model which is a proxy between the view that shows data, and the source model (TableModel), which stores this data.
    This proxy propagates data from the source to the view iff data meets the conditions stated in the filtering logic of
    this class.
    """

    def __init__(self, cols):
        super().__init__()
        # columns from which data will be taken into account for filtering
        self.cols = cols

    def filterAcceptsRow(self, sourceRow: int, sourceParent: QModelIndex):
        if sourceParent.isValid():
            source_model = self.sourceModel()

            # data will be concatenated from the data in the relevant columns into one string, which will be filtered.
            data = ''
            for col in self.cols:
                index = source_model.index(sourceRow, col, sourceParent)
                if col == 0:
                    filepath = str(source_model.data(index, Qt.UserRole))
                    data += get_show_name(filepath)
                    if 'error' in filepath:
                        data += '_pure'
                    if 'pseudo' in filepath:
                        data += '_pseudo'
                    if 'index' in filepath:
                        data += '_index'
                else:
                    data += str(source_model.data(index, Qt.DisplayRole))

            expression = self.filterRegularExpression()
            match = expression.match(data)
            return match.hasMatch()
        else:
            return True

    @staticmethod
    def filter(proxy, text):
        """
        Called upon user's filtering request, updates the regular expression that is used for automatic filtering of the
        model rows.
        """
        pattern = FilterModel.get_pattern(text)
        pattern = QRegularExpression(pattern, QRegularExpression.CaseInsensitiveOption)
        proxy.setFilterRegularExpression(pattern)

    @staticmethod
    def get_pattern(text):
        """
        User's request is divided on parts by '+', each part filters smth else like date, clustering purity or name.
        """
        splits = text.split('+')
        pattern = ''.join(['(?=.*' + split.strip() + ')' for split in splits])
        return pattern


class IconRenameDelegate(QStyledItemDelegate):
    """
    This class provides functionality of an icon +/- for expanding/shrinking rows as well as for renaming the
    rows. This delegate is supposed to work on the specific column. Higher hierarchy rows (groups of technologies)
    will be displaying an icon at that column, while on the next hierarchical level there will be the rows of data
    itself and the same col will hold the name of the dataset for which this delegate provides a renaming feature.
    """

    def __init__(self, ui, col, parent=None):
        super(IconRenameDelegate, self).__init__(parent)
        self._plus_icon = QIcon("icons\plus.png")
        self._minus_icon = QIcon("icons\minus.png")

        self.ui = ui
        self.col = col
        self.view = parent
        # to keep track and close the user input editors for renaming
        self.current_index = None

        self.view.doubleClicked.connect(self.cell_double_clicked)

    def initStyleOption(self, option, index):
        super(IconRenameDelegate, self).initStyleOption(option, index)
        if not index.parent().isValid():
            is_open = bool(option.state & QStyle.State_Open)
            option.features |= QStyleOptionViewItem.HasDecoration
            option.icon = self._minus_icon if is_open else self._plus_icon

        """
        this function is automatically invoked by a signal when the user submitted the new 
        desired name by double clicking and renaming cell in the table.
        """

    def setModelData(self, editor, model, index):
        if index.parent().isValid():
            new_name_suffix = editor.text()
            old_filepath = self.view.model().index(index.row(), self.col, index.parent()).data(Qt.UserRole)
            old_name_prefix = get_show_name(old_filepath)
            # nothing changed:
            if new_name_suffix == old_name_prefix:
                self.view.closePersistentEditor(self.current_index)
                return

            # empty name:
            if new_name_suffix == '':
                self.view.closePersistentEditor(self.current_index)
                self.ui.msg_box_with_title('New file name cannot be empty', 'Error!')
                return

            # check whether current file is being used:
            paths_in_work = [worker.input_path for worker in self.ui.cluster_workers + self.ui.error_workers]
            if old_filepath in paths_in_work:
                self.view.closePersistentEditor(self.current_index)
                self.ui.msg_box_with_title('This dataset is being used for a task, please wait until it finises to rename.', 'Error!')
                return

            # will check uniqueness as well:
            new_filepath = replace_suffix_path(old_filepath, new_name_suffix)
            if new_filepath == '':
                self.view.closePersistentEditor(self.current_index)
                self.ui.msg_box_with_title('Please choose a unique name for this dataset.', 'Error!')
                return

            # all is good, can rename:
            rename_json_entry(self.ui.recon_json, old_filepath, new_filepath, 'input_file')
            rename_json_entry(self.ui.clustering_res_json, old_filepath, new_filepath, 'name')
            # potential files for renaming:
            old_names = [old_filepath, get_shuffled_path_from_orig(old_filepath)]
            new_names = [new_filepath, get_shuffled_path_from_orig(new_filepath)]

            # pass over the clustering intermediate files:
            old_base_name = get_base_name(old_filepath)
            new_base_name = get_base_name(new_filepath)
            cluster_dir = Clustering.workdir()
            for file in files(cluster_dir):
                if old_base_name in file:  # found a candidate for renaming:
                    old_names.append(join_path(cluster_dir, file))
                    new_filename = join_path(cluster_dir, new_base_name + file.split(old_base_name)[1])
                    new_names.append(new_filename)

            for old, new in zip(old_names, new_names):
                rename_file(old, new)

            # rename filename in the tables:
            for table in self.ui.tables:
                model = TableModel.get_source_model(table)
                model.rename(old_filepath, new_filepath)

            # set up current input to a new filepath:
            self.ui.select_input(self.view, index)

            # close editor input line:
            self.view.closePersistentEditor(self.current_index)

        else:
            QStyledItemDelegate.setModelData(self, editor, model, index)

    def cell_double_clicked(self, index):
        if index.parent().isValid() and index.column() == self.col:
            self.view.openPersistentEditor(index)
            self.current_index = index


class ButtonDelegate(QStyledItemDelegate):
    """
    The cells in a column of a view for which this delegate will be used, will be acting as buttons.
    * paint() method is for drawing the buttons (those are not actually widgets, just pixmaps).
    * createEditor() is invoked upon user's mouse hovering event and will create the button widget for
    user interaction, the relevant slot for click event will be set in this function as well.
    """

    def __init__(self, ui, col, parent: QTreeView = None):
        super(ButtonDelegate, self).__init__(parent)
        self.ui = ui
        self.filepath_col = 0
        self.delegate_col = col
        self.button = QPushButton('', parent)
        self.view = parent
        self.current_index = None

        self.button.hide()
        self.view.entered.connect(self.cell_entered)
        self.view.setMouseTracking(True)

    def createEditor(self, parent, option, index):
        if index.parent().isValid():
            button = QPushButton(parent)
            self.set_button_view(button, option, index)
            filepath = self.view.model().index(index.row(), self.filepath_col, index.parent()).data(Qt.UserRole)
            button.clicked.connect(lambda: self.button_clicked(filepath))
            return button
        else:
            QStyledItemDelegate.createEditor(self, parent, option, index)

    def setEditorData(self, editor: QPushButton, index):
        if index.parent().isValid():
            editor.setProperty('data_value', index.data())
        else:
            QStyledItemDelegate.setEditorData(self, editor, index)

    def setModelData(self, editor, model, index):
        if index.parent().isValid():
            model.setData(index, editor.property('data_value'))
        else:
            QStyledItemDelegate.setModelData(self, editor, model, index)

    def paint(self, painter, option, index):
        if index.parent().isValid():
            self.set_button_view(self.button, option, index)
            pixmap = self.button.grab()
            painter.drawPixmap(option.rect.x(), option.rect.y(), pixmap)
        else:
            QStyledItemDelegate.paint(self, painter, option, index)

    def cell_entered(self, index):
        if self.current_index is not None:
            self.view.closePersistentEditor(self.current_index)
        if index.parent().isValid() and index.column() == self.delegate_col:
            self.view.openPersistentEditor(index)
            self.current_index = index

    def set_button_view(self, button, option, index):
        button.setText(index.data())
        button.setGeometry(option.rect)
        button.setFixedSize(option.rect.size())

    def button_clicked(self, data, info=None):
        pass


class DeleteButtonDelegate(ButtonDelegate):
    def __init__(self, ui, col, parent: QTreeView = None):
        super(DeleteButtonDelegate, self).__init__(ui, col, parent)

    def button_clicked(self, filepath, info=None):
        # stop active workers associated with file and return, files will be deleted with a termination signal afterwards:
        for worker in self.ui.error_workers + self.ui.cluster_workers + self.ui.recon_processes:
            if worker.input_path == filepath:
                worker.kill()
                return
        # auxiliary files that are connected to the original file being deleted.
        to_delete = [get_shuffled_path_from_orig(filepath)]
        # go through clustering files and delete relevant ones:
        base_path = get_base_name(filepath)
        cluster_dir = Clustering.workdir()
        for file in files(cluster_dir):
            if base_path in file:  # found an intermediate clustering dict file that needs to be deleted as well:
                to_delete.append(join_path(cluster_dir, file))

        # remove file from the table:
        TableModel.delete_file(self.ui.tables, filepath)
        # delete all the rest of associated files:
        for file in to_delete:
            TableModel.delete_file(None, file)
        # remove results from json files:
        update_json_file(self.ui.recon_json, None, {'input_file': filepath})
        update_json_file(self.ui.clustering_res_json, None, {'name': filepath})


class ShowButtonDelegate(ButtonDelegate):
    def __init__(self, ui, col, parent: QTreeView = None):
        super(ShowButtonDelegate, self).__init__(ui, col, parent)

    def createEditor(self, parent, option, index):
        if index.parent().isValid():
            button = QPushButton(parent)
            self.set_button_view(button, option, index)
            data = self.enable_button(button, index)
            button.clicked.connect(lambda: self.button_clicked(data))
            return button
        else:
            QStyledItemDelegate.createEditor(self, parent, option, index)

    def paint(self, painter, option, index):
        if index.parent().isValid():
            self.set_button_view(self.button, option, index)
            self.enable_button(self.button, index)
            pixmap = self.button.grab()
            painter.drawPixmap(option.rect.x(), option.rect.y(), pixmap)
        else:
            QStyledItemDelegate.paint(self, painter, option, index)

    def enable_button(self, button, index):
        pass

    @staticmethod
    def is_res_exist(results, data):
        exists = False
        with open(results) as rf:
            entries = json.load(rf)
        for entry in entries:
            for key in data:
                if entry[key] == data[key]:
                    exists = True
                else:
                    exists = False
                    break
            if exists:
                return True
        return False


class ClusterShowButtonDelegate(ShowButtonDelegate):
    def __init__(self, ui, col, parent: QTreeView = None):
        super(ClusterShowButtonDelegate, self).__init__(ui, col, parent)

    def enable_button(self, button, index):
        # get full path of the relevant file from the table model:
        filepath = self.view.model().index(index.row(), self.filepath_col, index.parent()).data(Qt.UserRole)
        # current state of the ui: affects which buttons to enable:
        info = {'name': filepath,
                'index': self.ui.clustering_index,
                'tech': self.ui.cluster_chosen_technology}
        enabled = ShowButtonDelegate.is_res_exist(self.ui.clustering_res_json, info)
        button.setEnabled(enabled)
        return info

    def button_clicked(self, data, info=None):
        found = False
        with open(self.ui.clustering_res_json) as rf:
            entries = json.load(rf)
        for entry in entries:
            for key in data:
                if entry[key] == data[key]:
                    found = True
                else:
                    found = False
                    break
            if found:  # show clustering results:
                self.ui.clustering_seconds_line.setText(str(entry['seconds']))
                self.ui.clustering_false_positive_line.setText(str(entry['false_positives']))
                self.ui.clustering_false_negative_line.setText(str(entry['false_negatives']))
                self.ui.clustering_strands_line.setText(str(entry['thrown_strands']))
                return


class ReconShowButtonDelegate(ShowButtonDelegate):
    def __init__(self, ui, col, parent: QTreeView = None):
        super(ReconShowButtonDelegate, self).__init__(ui, col, parent)

    def enable_button(self, button, index):
        # get full path of the relevant file from the table model:
        filepath = self.view.model().index(index.row(), self.filepath_col, index.parent()).data(Qt.UserRole)
        # current state of the ui: affects which buttons to enable:
        info = {'input_file': filepath,
                'algo': self.ui.recon_algo_dict[self.ui.reconstruction_algo]}
        enabled = ShowButtonDelegate.is_res_exist(self.ui.recon_json, info)
        button.setEnabled(enabled)
        return info

    def button_clicked(self, data, info=None):
        result_cols = ['input_file', 'output_file', 'error_rate', 'substitution_rate', 'deletion_rate',
                       'insertion_rate',
                       'success_rate', 'num_of_clusters', 'maj_test', 'hist_file']
        json_result = JsonReconResult()
        found = False
        with open(self.ui.recon_json) as rf:
            entries = json.load(rf)
        for entry in entries:
            for key in data:
                if entry[key] == data[key]:
                    found = True
                else:
                    found = False
                    break
            if found:  # show recon results:
                for col in result_cols:
                    json_result.__dict__[col] = entry[col]
                # create new window for displaying resutls:
                recon_result_window = ReconResult(json_result)
                self.ui.recon_result_windows.append(recon_result_window)
                recon_result_window.open_window()
                return


class StopButtonDelegate(ButtonDelegate):
    def __init__(self, ui, col, parent: QTreeView = None):
        super(StopButtonDelegate, self).__init__(ui, col, parent)

    def createEditor(self, parent, option, index):
        button = QPushButton(parent)
        self.set_button_view(button, option, index)
        self.enable_button(button, index)
        filepath = self.view.model().index(index.row(), self.filepath_col, index.parent()).data(Qt.UserRole)
        info = {}
        if self.view == self.ui.clustering_progress_table:
            # tech + index:
            info['cluster_tech'] = self.view.model().index(index.row(), self.ui.cluster_progress_cols['tech'], index.parent()).data()
            info['cluster_index'] = self.view.model().index(index.row(), self.ui.cluster_progress_cols['index'], index.parent()).data()
        elif self.view == self.ui.recon_progress_table:
            # algo:
            info['algo'] = self.view.model().index(index.row(), self.ui.recon_progress_cols['algo'], index.parent()).data()
        button.clicked.connect(lambda: self.button_clicked(filepath, info))
        return button

    def paint(self, painter, option, index):
        self.set_button_view(self.button, option, index)
        self.enable_button(self.button, index)
        pixmap = self.button.grab()
        painter.drawPixmap(option.rect.x(), option.rect.y(), pixmap)

    def enable_button(self, button, index):
        enabled_str = self.view.model().index(index.row(), self.delegate_col, index.parent()).data(Qt.UserRole)
        enabled = True if enabled_str == 'enabled' else False
        button.setEnabled(enabled)

    def button_clicked(self, filepath, info=None):
        workers = []
        if self.view == self.ui.clustering_progress_table:
            workers = self.ui.cluster_workers
        elif self.view == self.ui.recon_progress_table:
            workers = self.ui.recon_processes
        # stop worker:
        for worker in workers:
            if worker.input_path == filepath:
                to_kill = True
                for key in info:
                    to_kill = to_kill * str(getattr(worker, key)) == info[key]
                if to_kill:
                    worker.kill()

    def cell_entered(self, index):
        if index.column() == self.delegate_col:
            if self.current_index is not None:
                self.view.closePersistentEditor(self.current_index)
            self.view.openPersistentEditor(index)
            self.current_index = index


class ProgressBarDelegate(QStyledItemDelegate):
    """
    Cells of column for which this delegate will be used, will act as a progressbars.
    * paint() will repaint the relevant cell upon update of the integer value in the cell which acts as
    width. This value is updated with a signal from process or worker.
    """

    def __init__(self, parent=None, expandable=True):
        super(ProgressBarDelegate, self).__init__(parent)
        self.view = parent
        self.expandable = expandable

    def paint(self, painter, option, index):
        if self.expandable and index.parent().isValid() or not self.expandable:
            value = int(index.data())
            rect = QRect(option.rect)
            rect.setWidth(int(value * rect.width() / 100))
            painter.fillRect(rect, QColor('lime').darker(value + 60))
            painter.drawText(option.rect, Qt.AlignmentFlag.AlignCenter, index.data())


class TableModel(QStandardItemModel):
    """
    Main building fundament of representing info as a table. This is not a table itself that is seen in the app,
    what we see is a view, this class is a model that holds the data in a hierarchical order for the view
    for display.
    """

    def __init__(self, horizontalLabels, view, parent=None):
        super(TableModel, self).__init__(parent)
        self.view = view
        self.setColumnCount(len(horizontalLabels))
        self.setHorizontalHeaderLabels(horizontalLabels)

    def create_item(self, data=None):
        if data:
            item = QStandardItem(data)
        else:
            item = QStandardItem()
        item.setEditable(False)
        return item

    def add_group(self, group_name):
        item = self.create_item(group_name)
        root = self.invisibleRootItem()
        rows = root.rowCount()
        for col, it in enumerate([item]):
            root.setChild(rows, col, it)

    def get_group_from_filepath(self, filepath):
        group_name = get_tech_from_path(filepath)
        group = self.get_group_from_group_name(group_name)
        return group

    def get_group_from_group_name(self, group_name):
        root = self.invisibleRootItem()
        rows = root.rowCount()
        for row in range(rows):
            if self.item(row, 0).text() == group_name:
                return self.item(row, 0)

    def append_row_to_group(self, row_elements, secondary_info=None, toExpand=False):
        """
        secondary_info is provided in case there is need to check uniqueness of rows in the table,
        also it is an indication to the fact that row is going to be added to the root directly (not to a group).
        In other words, secondary_info means, append to progress table which is not hierarchical.
        The returned value indicates rather a row was really appended.
        """
        col_items = [QStandardItem() for _ in range(len(row_elements))]
        filepath = row_elements[0]['hidden']

        if secondary_info is None:
            group = self.get_group_from_filepath(get_base_name(filepath))
            row = None
        else:
            group = self.invisibleRootItem()
            row = self.get_row_in_group(group, filepath, secondary_info)
        if row is None:  # no such row in the table, lets add:
            for col in row_elements:
                data = row_elements[col]
                item = self.create_item(data['data'])
                if 'hidden' in data:
                    item.setData(data['hidden'], Qt.UserRole)
                col_items[col] = item
            group.appendRow(col_items)
            # expand the relevant group for the user to highlight where the row was added:
            self.view.setExpanded(group.index(), toExpand)
            return True
        return False

    def remove_row(self, filepath, in_root=False, secondary_info=None):
        group = self.invisibleRootItem() if in_root else self.get_group_from_filepath(filepath)
        row = self.get_row_in_group(group, filepath, secondary_info)
        if row is not None:
            self.view.selectionModel().clear()
            self.removeRow(row, group.index())

    def get_row_in_group(self, group_item, row_name, secondary_info=None):
        rows = group_item.rowCount()
        for row in range(rows):
            if group_item.child(row, 0).data(Qt.UserRole) == row_name:
                if secondary_info:
                    found = True
                    for data, col in zip(secondary_info['data'], secondary_info['cols']):
                        found = found and group_item.child(row, col).text() == str(data)
                    if found:
                        return row
                else:
                    return row
        return None

    def update_cell(self, value, filepath, col, in_root=False, secondary_info=None, role=Qt.DisplayRole):
        """
        updates specified cell in an existing row in the table.
        """
        group = self.invisibleRootItem() if in_root else self.get_group_from_filepath(filepath)
        row = self.get_row_in_group(group, filepath, secondary_info)

        item = group.child(row, col)
        item.setData(str(value), role)
        if row is not None:
            group.setChild(row, col, item)

    def rename(self, old_filepath, new_filepath, in_root=False):
        group = self.invisibleRootItem() if in_root else self.get_group_from_filepath(old_filepath)
        rows = group.rowCount()
        role = Qt.UserRole
        for row in range(rows):
            item = group.child(row, 0)
            filepath = item.data(role=role)
            if filepath == old_filepath:
                item.setData(new_filepath, role)
                item.setText(get_show_name(new_filepath))
                group.setChild(row, 0, item)

    @pyqtSlot(QModelIndex)
    def on_clicked(self, index, view, expandable):
        if expandable and not index.parent().isValid() and index.column() == 0:
            view.setExpanded(index, not view.isExpanded(index))

    @staticmethod
    def get_source_model(table):
        """In case there is a proxy model, all the relevant info is stored in most lower model layer"""
        model = table.model()
        # remove all layers of models to get to the initial source:
        while hasattr(model, 'sourceModel'):
            model = model.sourceModel()
        return model

    @staticmethod
    def set_proxy(table, model, cols):
        cols = [cols['name'], cols['time']]
        proxy = FilterModel(cols)
        proxy.setSourceModel(model)
        table.setModel(proxy)
        return proxy

    @staticmethod
    def delete_file(tables, filepath):
        # delete rows associated with the file from all the tables.
        if tables is not None:
            for table in tables:
                model = TableModel.get_source_model(table)
                model.remove_row(filepath)
        # delete the file itself:
        delete(filepath)
